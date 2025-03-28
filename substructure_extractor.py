#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
内部自由度に基づく部分構造自動抽出アルゴリズムを提供するモジュール
"""

import os

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def identify_rigid_substructure(mol_file, target_atom_idx):
    """
    置換対象原子を中心とした内部自由度の少ない部分構造を自動的に特定する
    
    Parameters
    ----------
    mol_file : str
        分子のSDFファイルパス
    target_atom_idx : int
        置換対象原子のインデックス
    
    Returns
    -------
    list
        抽出すべき原子のインデックスリスト
    """
    # 分子の読み込み
    mol_supplier = Chem.SDMolSupplier(mol_file)
    mol = mol_supplier[0]
    
    if mol is None:
        raise ValueError(f"分子ファイル {mol_file} を読み込めませんでした")
    
    # 初期化：置換対象原子を含める
    selected_atoms = set([target_atom_idx])
    
    # 環構造の検出と保持
    ring_info = mol.GetRingInfo()
    for ring_atoms in ring_info.AtomRings():
        if target_atom_idx in ring_atoms:
            selected_atoms.update(ring_atoms)
    
    # 共役系の検出と保持
    conjugated_atoms = identify_conjugated_system(mol, selected_atoms)
    selected_atoms.update(conjugated_atoms)
    
    # 回転自由度の少ない結合を持つ原子を追加
    rigid_neighbors = identify_rigid_bonds(mol, selected_atoms)
    selected_atoms.update(rigid_neighbors)
    
    # 内部自由度を最小化する境界を決定
    optimal_boundary = optimize_boundary(mol, selected_atoms)
    
    return list(optimal_boundary)


def identify_conjugated_system(mol, seed_atoms):
    """
    指定した原子群から連続する共役系を特定する
    
    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        分子
    seed_atoms : set
        開始点となる原子のインデックスセット
    
    Returns
    -------
    set
        共役系に含まれる原子のインデックスセット
    """
    conjugated_atoms = set()
    
    # 共役系の特定ロジック
    # - 芳香環の検出
    # - 二重結合、三重結合の検出
    # - sp2炭素の連鎖の検出
    
    for atom_idx in seed_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        
        # 芳香族原子の場合
        if atom.GetIsAromatic():
            # 同じ芳香環に属する原子を追加
            for bond in atom.GetBonds():
                other_atom = bond.GetOtherAtom(atom)
                if other_atom.GetIsAromatic():
                    conjugated_atoms.add(other_atom.GetIdx())
        
        # sp2またはsp混成の原子の場合
        hybridization = atom.GetHybridization()
        if hybridization == 2 or hybridization == 1:  # SP2=2, SP=1
            # 結合している同様の混成状態の原子を追加
            for bond in atom.GetBonds():
                bond_type = bond.GetBondType()
                if bond_type == 2 or bond_type == 3:  # DOUBLE=2, TRIPLE=3
                    other_atom = bond.GetOtherAtom(atom)
                    conjugated_atoms.add(other_atom.GetIdx())
    
    return conjugated_atoms


def identify_rigid_bonds(mol, seed_atoms):
    """
    回転自由度の少ない結合を持つ原子を特定する
    
    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        分子
    seed_atoms : set
        開始点となる原子のインデックスセット
    
    Returns
    -------
    set
        回転自由度の少ない結合を持つ原子のインデックスセット
    """
    rigid_neighbors = set()
    
    for atom_idx in seed_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        
        for bond in atom.GetBonds():
            # 二重結合、三重結合、または環内結合の場合
            bond_type = bond.GetBondType()
            if (bond_type == 2 or  # DOUBLE=2
                bond_type == 3 or  # TRIPLE=3
                bond.IsInRing()):
                other_atom = bond.GetOtherAtom(atom)
                rigid_neighbors.add(other_atom.GetIdx())
    
    return rigid_neighbors


def optimize_boundary(mol, selected_atoms):
    """
    内部自由度を最小化する境界を決定する
    
    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        分子
    selected_atoms : set
        選択された原子のインデックスセット
    
    Returns
    -------
    set
        最適化された原子のインデックスセット
    """
    # 現在の選択に基づく内部自由度を計算
    current_dof = calculate_internal_dof(mol, selected_atoms)
    
    # 境界原子（選択された原子に直接結合しているが、選択されていない原子）を特定
    boundary_atoms = set()
    for atom_idx in selected_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in selected_atoms:
                boundary_atoms.add(neighbor_idx)
    
    # 各境界原子を追加した場合の内部自由度の変化を評価
    optimal_selection = set(selected_atoms)
    
    for boundary_idx in boundary_atoms:
        test_selection = selected_atoms.union([boundary_idx])
        test_dof = calculate_internal_dof(mol, test_selection)
        
        # 内部自由度の増加が最小限であれば追加
        if test_dof - current_dof <= 1:  # 閾値は調整可能
            optimal_selection.add(boundary_idx)
            current_dof = test_dof
    
    return optimal_selection


def calculate_internal_dof(mol, atom_indices):
    """
    原子群の内部自由度を計算する
    
    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        分子
    atom_indices : set
        原子のインデックスセット
    
    Returns
    -------
    int
        内部自由度
    """
    # 原子数
    n_atoms = len(atom_indices)
    
    # 基本的な内部自由度: 3N-6（N=原子数）
    basic_dof = 3 * n_atoms - 6 if n_atoms >= 3 else 0
    
    # 環構造による拘束
    ring_constraints = 0
    ring_info = mol.GetRingInfo()
    for ring_atoms in ring_info.AtomRings():
        ring_atoms_set = set(ring_atoms)
        if ring_atoms_set.issubset(atom_indices):
            # 環ごとに自由度が減少
            if len(ring_atoms_set) > 3:  # 3員環以上
                ring_constraints += len(ring_atoms_set) - 3
    
    # 二重結合、三重結合による拘束
    bond_constraints = 0
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        if begin_idx in atom_indices and end_idx in atom_indices:
            bond_type = bond.GetBondType()
            if bond_type == 2:  # DOUBLE=2
                bond_constraints += 1
            elif bond_type == 3:  # TRIPLE=3
                bond_constraints += 2
    
    # 実効的な内部自由度
    effective_dof = basic_dof - ring_constraints - bond_constraints
    
    return max(0, effective_dof)


def extract_substructure(mol_file, atom_indices):
    """
    指定した原子インデックスの集合を含む部分構造を抽出する
    
    Parameters
    ----------
    mol_file : str
        分子のSDFファイルパス
    atom_indices : list
        抽出する原子のインデックスリスト
    
    Returns
    -------
    rdkit.Chem.rdchem.Mol
        抽出された部分構造
    """
    # 分子の読み込み
    mol_supplier = Chem.SDMolSupplier(mol_file)
    mol = mol_supplier[0]
    
    if mol is None:
        raise ValueError(f"分子ファイル {mol_file} を読み込めませんでした")
    
    # 部分構造を抽出するためにEditable Molを使用
    emol = Chem.EditableMol(Chem.Mol())
    
    # 原子マッピングを作成
    atom_mapping = {}
    
    # 選択した原子を追加
    for old_idx in atom_indices:
        atom = mol.GetAtomWithIdx(old_idx)
        # 原子番号を取得
        atomic_num = atom.GetAtomicNum()
        # 新しい原子を追加
        new_idx = emol.AddAtom(Chem.Atom(atomic_num))
        atom_mapping[old_idx] = new_idx
    
    # 結合を追加
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        if begin_idx in atom_indices and end_idx in atom_indices:
            emol.AddBond(
                atom_mapping[begin_idx],
                atom_mapping[end_idx],
                bond.GetBondType()
            )
    
    # Editable Molから通常の分子に変換
    substructure = emol.GetMol()
    
    # 3D座標を設定
    if mol.GetNumConformers() > 0:
        conf = Chem.Conformer(substructure.GetNumAtoms())
        for old_idx, new_idx in atom_mapping.items():
            old_pos = mol.GetConformer().GetAtomPosition(old_idx)
            conf.SetAtomPosition(new_idx, old_pos)
        
        substructure.AddConformer(conf)
    
    # 切断された結合の末端に水素を追加して安定化
    substructure = Chem.AddHs(substructure)
    
    # 3D座標を最適化
    try:
        AllChem.EmbedMolecule(substructure)
        AllChem.UFFOptimizeMolecule(substructure)
    except Exception as e:
        print(f"警告: 3D座標の最適化に失敗しました: {e}")
    
    return substructure


def replace_atom(mol, atom_idx, new_atom_symbol):
    """
    指定した原子を別の原子に置換する
    
    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        分子
    atom_idx : int
        置換する原子のインデックス
    new_atom_symbol : str
        新しい原子のシンボル（例: 'N', 'O'）
    
    Returns
    -------
    rdkit.Chem.rdchem.Mol
        置換後の分子
    """
    # 新しい分子を作成
    emol = Chem.EditableMol(Chem.Mol())
    
    # 原子マッピングを作成
    atom_mapping = {}
    
    # 原子を追加
    for idx in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(idx)
        
        # 置換対象の原子の場合
        if idx == atom_idx:
            # 新しい原子番号を取得
            if new_atom_symbol == 'C': atomic_num = 6
            elif new_atom_symbol == 'N': atomic_num = 7
            elif new_atom_symbol == 'O': atomic_num = 8
            elif new_atom_symbol == 'F': atomic_num = 9
            elif new_atom_symbol == 'P': atomic_num = 15
            elif new_atom_symbol == 'S': atomic_num = 16
            elif new_atom_symbol == 'Cl': atomic_num = 17
            elif new_atom_symbol == 'Br': atomic_num = 35
            elif new_atom_symbol == 'I': atomic_num = 53
            else:
                raise ValueError(f"未対応の元素記号: {new_atom_symbol}")
            
            # 新しい原子を追加
            new_idx = emol.AddAtom(Chem.Atom(atomic_num))
        else:
            # 元の原子をそのまま追加
            new_idx = emol.AddAtom(Chem.Atom(atom.GetAtomicNum()))
        
        atom_mapping[idx] = new_idx
    
    # 結合を追加
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        emol.AddBond(
            atom_mapping[begin_idx],
            atom_mapping[end_idx],
            bond.GetBondType()
        )
    
    # Editable Molから通常の分子に変換
    replaced_mol = emol.GetMol()
    
    # 3D座標を設定
    if mol.GetNumConformers() > 0:
        conf = Chem.Conformer(replaced_mol.GetNumAtoms())
        for old_idx, new_idx in atom_mapping.items():
            old_pos = mol.GetConformer().GetAtomPosition(old_idx)
            conf.SetAtomPosition(new_idx, old_pos)
        
        replaced_mol.AddConformer(conf)
    
    # 形式電荷や芳香族性などの特性を設定
    for old_idx, new_idx in atom_mapping.items():
        old_atom = mol.GetAtomWithIdx(old_idx)
        new_atom = replaced_mol.GetAtomWithIdx(new_idx)
        
        # 形式電荷を設定
        new_atom.SetFormalCharge(old_atom.GetFormalCharge())
        
        # 芳香族性を設定
        new_atom.SetIsAromatic(old_atom.GetIsAromatic())
        
        # キラリティを設定
        new_atom.SetChiralTag(old_atom.GetChiralTag())
    
    # 暗黙的な原子価を計算
    Chem.SanitizeMol(replaced_mol)
    
    return replaced_mol


def process_molecule_for_msmd(mol_file, target_atom_idx, new_atom_symbol):
    """
    分子の指定した原子を置換し、部分構造を抽出して量子化学計算用のファイルを作成する
    
    Parameters
    ----------
    mol_file : str
        分子のSDFファイルパス
    target_atom_idx : int
        置換対象原子のインデックス
    new_atom_symbol : str
        新しい原子のシンボル（例: 'N', 'O'）
    
    Returns
    -------
    tuple
        (元の部分構造のファイルパス, 置換後の部分構造のファイルパス)
    """
    # 出力ディレクトリを作成
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # 内部自由度に基づく部分構造の自動抽出
    atom_indices = identify_rigid_substructure(mol_file, target_atom_idx)
    
    # 部分構造の抽出
    substructure = extract_substructure(mol_file, atom_indices)
    
    # 元の部分構造を保存
    original_file = os.path.join(output_dir, "original_substructure.sdf")
    writer = Chem.SDWriter(original_file)
    writer.write(substructure)
    writer.close()
    
    # 置換対象原子のインデックスを新しい部分構造内のインデックスに変換
    # （この例では単純化のため、元のインデックスと同じと仮定）
    new_target_idx = atom_indices.index(target_atom_idx)
    
    # 原子を置換
    replaced_substructure = replace_atom(substructure, new_target_idx, new_atom_symbol)
    
    # 置換後の部分構造を保存
    replaced_file = os.path.join(output_dir, "replaced_substructure.sdf")
    writer = Chem.SDWriter(replaced_file)
    writer.write(replaced_substructure)
    writer.close()
    
    print(f"元の部分構造を {original_file} に保存しました")
    print(f"置換後の部分構造を {replaced_file} に保存しました")
    
    return original_file, replaced_file


def main():
    """
    メイン関数
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='内部自由度に基づく部分構造自動抽出')
    parser.add_argument('--mol', required=True, help='分子のSDFファイルパス')
    parser.add_argument('--atom-idx', type=int, required=True, help='置換対象原子のインデックス')
    parser.add_argument('--new-atom', required=True, help='新しい原子のシンボル（例: N, O）')
    args = parser.parse_args()
    
    # 分子の処理
    original_file, replaced_file = process_molecule_for_msmd(
        args.mol, args.atom_idx, args.new_atom
    )
    
    print("次のステップ: 生成された部分構造に対してquantitative inverse MSMDを実行してください")


if __name__ == "__main__":
    main()