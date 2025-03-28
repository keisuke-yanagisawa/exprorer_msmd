#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
MCS計算と非共通部分抽出機能を提供するモジュール
"""

import os

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS


def extract_non_common_parts(mol1_file, mol2_file):
    """
    2つの化合物間のMCSを計算し、非共通部分を抽出する

    Parameters
    ----------
    mol1_file : str
        1つ目の化合物のSDFファイルパス
    mol2_file : str
        2つ目の化合物のSDFファイルパス

    Returns
    -------
    dict
        各化合物の非共通部分の情報
        {
            'mol1': {'mol': rdkit.Chem.rdchem.Mol, 'atom_indices': list},
            'mol2': {'mol': rdkit.Chem.rdchem.Mol, 'atom_indices': list}
        }
    """
    # 分子の読み込み
    mol1 = Chem.SDMolSupplier(mol1_file)[0]
    mol2 = Chem.SDMolSupplier(mol2_file)[0]
    
    if mol1 is None:
        raise ValueError(f"分子ファイル {mol1_file} を読み込めませんでした")
    if mol2 is None:
        raise ValueError(f"分子ファイル {mol2_file} を読み込めませんでした")
    
    # MCSの計算
    mcs_result = rdFMCS.FindMCS(
        [mol1, mol2],
        completeRingsOnly=True,
        ringMatchesRingOnly=True,
        matchValences=True
    )
    
    # MCSのSMARTSパターンを取得
    mcs_smarts = mcs_result.smartsString
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)
    
    # 各分子でMCSに対応する原子インデックスを取得
    match1 = mol1.GetSubstructMatch(mcs_mol)
    match2 = mol2.GetSubstructMatch(mcs_mol)
    
    # 非共通部分の原子インデックスを取得
    all_atoms1 = set(range(mol1.GetNumAtoms()))
    all_atoms2 = set(range(mol2.GetNumAtoms()))
    
    non_common_atoms1 = list(all_atoms1 - set(match1))
    non_common_atoms2 = list(all_atoms2 - set(match2))
    
    # 非共通部分を含むサブ構造を抽出
    non_common_mol1 = get_substructure(mol1, non_common_atoms1, match1)
    non_common_mol2 = get_substructure(mol2, non_common_atoms2, match2)
    
    return {
        'mol1': {'mol': non_common_mol1, 'atom_indices': non_common_atoms1},
        'mol2': {'mol': non_common_mol2, 'atom_indices': non_common_atoms2}
    }


def get_substructure(mol, atom_indices, mcs_indices):
    """
    指定された原子インデックスに基づいてサブ構造を抽出
    MCSとの接続点も含める

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        元の分子
    atom_indices : list
        抽出する原子のインデックスリスト
    mcs_indices : list
        MCSに対応する原子のインデックスリスト

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        抽出されたサブ構造
    """
    # 非共通部分の原子に加えて、MCSとの接続点も含める
    connection_atoms = []
    
    for atom_idx in atom_indices:
        atom = mol.GetAtomWithIdx(atom_idx)
        
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            
            if neighbor_idx in mcs_indices and neighbor_idx not in connection_atoms:
                connection_atoms.append(neighbor_idx)
    
    # 抽出する原子インデックスを結合
    extraction_indices = atom_indices + connection_atoms
    
    # サブ構造を抽出
    # RDKitでは直接サブ構造を抽出する方法が限られているため、
    # 元の分子をコピーして、抽出する原子以外を削除する方法を使用
    
    # 元の分子をコピー
    substructure = Chem.Mol(mol)
    
    # 抽出する原子以外を削除するためのマップを作成
    atoms_to_keep = {idx: idx for idx in extraction_indices}
    
    # 新しい分子を作成
    edit_mol = Chem.EditableMol(Chem.Mol())
    
    # 抽出する原子を追加
    new_idx_map = {}
    for old_idx in extraction_indices:
        atom = mol.GetAtomWithIdx(old_idx)
        new_idx = edit_mol.AddAtom(atom)
        new_idx_map[old_idx] = new_idx
    
    # 結合を追加
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        if begin_idx in extraction_indices and end_idx in extraction_indices:
            edit_mol.AddBond(
                new_idx_map[begin_idx],
                new_idx_map[end_idx],
                bond.GetBondType()
            )
    
    # 編集可能な分子から通常の分子に変換
    substructure = edit_mol.GetMol()
    
    # 3D座標を設定
    conf = Chem.Conformer(len(extraction_indices))
    for old_idx, new_idx in new_idx_map.items():
        old_pos = mol.GetConformer().GetAtomPosition(old_idx)
        conf.SetAtomPosition(new_idx, old_pos)
    
    substructure.AddConformer(conf)
    
    return substructure


def visualize_mcs_and_non_common_parts(mol1_file, mol2_file, output_dir="."):
    """
    MCSと非共通部分を可視化する

    Parameters
    ----------
    mol1_file : str
        1つ目の化合物のSDFファイルパス
    mol2_file : str
        2つ目の化合物のSDFファイルパス
    output_dir : str, optional
        出力ディレクトリ, by default "."
    """
    # 分子の読み込み
    mol1 = Chem.SDMolSupplier(mol1_file)[0]
    mol2 = Chem.SDMolSupplier(mol2_file)[0]
    
    # MCSの計算
    mcs_result = rdFMCS.FindMCS(
        [mol1, mol2],
        completeRingsOnly=True,
        ringMatchesRingOnly=True,
        matchValences=True
    )
    
    # MCSのSMARTSパターンを取得
    mcs_smarts = mcs_result.smartsString
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)
    
    # 各分子でMCSに対応する原子インデックスを取得
    match1 = mol1.GetSubstructMatch(mcs_mol)
    match2 = mol2.GetSubstructMatch(mcs_mol)
    
    # MCSを強調表示
    mol1_highlight = Chem.Mol(mol1)
    mol2_highlight = Chem.Mol(mol2)
    
    # 非共通部分を抽出
    non_common_parts = extract_non_common_parts(mol1_file, mol2_file)
    
    # 結果を保存
    os.makedirs(output_dir, exist_ok=True)
    
    # 非共通部分を保存
    writer1 = Chem.SDWriter(os.path.join(output_dir, "non_common_mol1.sdf"))
    writer1.write(non_common_parts['mol1']['mol'])
    writer1.close()
    
    writer2 = Chem.SDWriter(os.path.join(output_dir, "non_common_mol2.sdf"))
    writer2.write(non_common_parts['mol2']['mol'])
    writer2.close()
    
    print(f"非共通部分を {output_dir} に保存しました")
    print(f"  mol1の非共通部分: non_common_mol1.sdf")
    print(f"  mol2の非共通部分: non_common_mol2.sdf")


def main():
    """
    メイン関数
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='MCS計算と非共通部分抽出ツール')
    parser.add_argument('--mol1', required=True, help='1つ目の化合物のSDFファイルパス')
    parser.add_argument('--mol2', required=True, help='2つ目の化合物のSDFファイルパス')
    parser.add_argument('--output-dir', default=".", help='出力ディレクトリ')
    args = parser.parse_args()
    
    # MCSと非共通部分を可視化
    visualize_mcs_and_non_common_parts(args.mol1, args.mol2, args.output_dir)
    
    # 非共通部分を抽出
    non_common_parts = extract_non_common_parts(args.mol1, args.mol2)
    
    print(f"mol1の非共通部分の原子数: {len(non_common_parts['mol1']['atom_indices'])}")
    print(f"mol2の非共通部分の原子数: {len(non_common_parts['mol2']['atom_indices'])}")


if __name__ == "__main__":
    main()