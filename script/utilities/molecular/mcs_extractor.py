#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
MCS計算と非共通部分抽出機能を提供するモジュール
"""

import os

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS


def create_substructure(mol, extraction_indices):
    """
    指定された原子インデックスに基づいてサブ構造を作成

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        元の分子
    extraction_indices : list
        抽出する原子のインデックスリスト

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        抽出されたサブ構造
    """
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


def extract_non_common_parts(mol1_file, mol2_file, only_non_common_rings=False, separate_rings=False):
    """
    2つの化合物間のMCSを計算し、非共通部分を抽出する

    Parameters
    ----------
    mol1_file : str
        1つ目の化合物のSDFファイルパス
    mol2_file : str
        2つ目の化合物のSDFファイルパス
    only_non_common_rings : bool, optional
        Trueの場合、複合環のうち非共通な環のみをプローブ化する, by default False
    separate_rings : bool, optional
        Trueの場合、非共通な環を別々のプローブとして抽出する, by default False

    Returns
    -------
    dict
        各化合物の非共通部分の情報
        separate_rings=Falseの場合:
        {
            'mol1': {'mol': rdkit.Chem.rdchem.Mol, 'atom_indices': list},
            'mol2': {'mol': rdkit.Chem.rdchem.Mol, 'atom_indices': list}
        }
        separate_rings=Trueの場合:
        {
            'mol1': [{'mol': rdkit.Chem.rdchem.Mol, 'atom_indices': list}, ...],
            'mol2': [{'mol': rdkit.Chem.rdchem.Mol, 'atom_indices': list}, ...]
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
    non_common_mol1 = get_substructure(mol1, non_common_atoms1, match1, only_non_common_rings, separate_rings)
    non_common_mol2 = get_substructure(mol2, non_common_atoms2, match2, only_non_common_rings, separate_rings)
    
    # separate_ringsがTrueの場合、複数のサブ構造を返す
    if separate_rings and only_non_common_rings:
        # 既に環ごとに分割されたサブ構造のリストが返されている
        return {
            'mol1': non_common_mol1,  # リスト
            'mol2': non_common_mol2   # リスト
        }
    else:
        # 単一のサブ構造が返されている
        return {
            'mol1': {'mol': non_common_mol1, 'atom_indices': non_common_atoms1},
            'mol2': {'mol': non_common_mol2, 'atom_indices': non_common_atoms2}
        }


def get_substructure(mol, atom_indices, mcs_indices, only_non_common_rings=False, separate_rings=False):
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
    only_non_common_rings : bool, optional
        Trueの場合、複合環のうち非共通な環のみを抽出する, by default False
    separate_rings : bool, optional
        Trueの場合、非共通な環を別々のプローブとして抽出する, by default False

    Returns
    -------
    rdkit.Chem.rdchem.Mol or list
        separate_rings=Falseの場合: 抽出されたサブ構造
        separate_rings=Trueの場合: 抽出されたサブ構造のリスト
    """
    # 環構造の情報を取得
    # 対称性を考慮した環のセットを取得
    try:
        from rdkit.Chem.rdmolops import GetSymmSSSR
        ring_atoms = GetSymmSSSR(mol)
    except:
        # GetSymmSSSRが利用できない場合は、通常の環情報を使用
        ring_info = mol.GetRingInfo()
        ring_atoms = ring_info.AtomRings()
    
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
    
    # 環ごとに分割するかどうかのフラグ
    if separate_rings and only_non_common_rings:
        # 非共通部分に含まれる環を特定
        non_common_rings = []
        common_rings = []
        
        # 環を共通部分と非共通部分に分類
        for ring in ring_atoms:
            # 環の原子が非共通部分に含まれるかチェック
            if any(atom_idx in atom_indices for atom_idx in ring):
                # 環全体が共通部分に含まれる場合は共通環とする
                if all(atom_idx in mcs_indices for atom_idx in ring):
                    common_rings.append(ring)
                # それ以外は非共通環とする
                else:
                    non_common_rings.append(ring)
        
        # 各環ごとに別々のサブ構造を作成
        if non_common_rings:
            # 各環ごとのサブ構造のリストを作成
            substructures = []
            
            for ring in non_common_rings:
                # 環の原子のうち、非共通部分に含まれる原子のみを追加
                ring_atom_indices = set()
                for atom_idx in ring:
                    if atom_idx in atom_indices:
                        ring_atom_indices.add(atom_idx)
                
                # 接続点を計算（環の原子に隣接する共通部分の原子のみを含める）
                ring_connection_atoms = []
                for atom_idx in ring_atom_indices:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        neighbor_idx = neighbor.GetIdx()
                        if neighbor_idx in mcs_indices and neighbor_idx not in ring_connection_atoms:
                            ring_connection_atoms.append(neighbor_idx)
                
                # 環の原子と接続点を抽出
                if ring_atom_indices:
                    ring_extraction_indices = list(ring_atom_indices) + ring_connection_atoms
                    
                    # サブ構造を作成
                    substructure = create_substructure(mol, ring_extraction_indices)
                    
                    # 環のSMILES文字列を取得
                    try:
                        ring_smiles = Chem.MolToSmiles(substructure)
                        
                        # 環のSMILES文字列が長すぎる場合（複数の環を含む場合）は、
                        # 単純な環のみを抽出する
                        if len(ring_smiles) > 30 or ring_smiles.count('c') > 10 or 'c2c' in ring_smiles or 'c3c' in ring_smiles:
                            print(f"警告: 環が複雑すぎます: {ring_smiles}")
                            continue
                        
                        # 環のSMILES文字列に複数の環を示す特徴がある場合は、
                        # 単純な環のみを抽出する
                        if ring_smiles.count('1') > 2 or ring_smiles.count('2') > 2 or ring_smiles.count('3') > 2:
                            print(f"警告: 環が複数含まれています: {ring_smiles}")
                            continue
                    except:
                        print("警告: 環のSMILES文字列を取得できませんでした")
                        continue
                    
                    substructures.append({
                        'mol': substructure,
                        'atom_indices': list(ring_atom_indices)
                    })
            
            # 複数のサブ構造を返す
            return substructures
            
        else:
            # 非共通な環が見つからない場合は、元の非共通部分を使用
            substructure = create_substructure(mol, extraction_indices)
            return {'mol': substructure, 'atom_indices': atom_indices}
    
    # only_non_common_ringsがTrueの場合、複合環のうち非共通な環のみを抽出
    elif only_non_common_rings:
        # 非共通部分に含まれる環を特定
        non_common_rings = []
        common_rings = []
        
        # 環を共通部分と非共通部分に分類
        for ring in ring_atoms:
            # 環の原子が非共通部分に含まれるかチェック
            if any(atom_idx in atom_indices for atom_idx in ring):
                # 環全体が共通部分に含まれる場合は共通環とする
                if all(atom_idx in mcs_indices for atom_idx in ring):
                    common_rings.append(ring)
                # それ以外は非共通環とする
                else:
                    non_common_rings.append(ring)
        
        # 非共通部分の環に含まれる原子のみを抽出
        if non_common_rings:
            # 非共通部分の環に含まれる原子を集める
            ring_atom_indices = set()
            for ring in non_common_rings:
                # 環の原子のうち、非共通部分に含まれる原子のみを追加
                for atom_idx in ring:
                    if atom_idx in atom_indices:
                        ring_atom_indices.add(atom_idx)
            
            # 非共通部分の環に含まれる原子と接続点を抽出
            if ring_atom_indices:
                # 接続点を再計算（非共通環の原子に隣接する共通部分の原子のみを含める）
                new_connection_atoms = []
                for atom_idx in ring_atom_indices:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        neighbor_idx = neighbor.GetIdx()
                        if neighbor_idx in mcs_indices and neighbor_idx not in new_connection_atoms:
                            new_connection_atoms.append(neighbor_idx)
                
                extraction_indices = list(ring_atom_indices) + new_connection_atoms
            else:
                # 非共通な環が見つからない場合は、元の非共通部分を使用
                pass
    
    # サブ構造を作成して返す
    substructure = create_substructure(mol, extraction_indices)
    return {'mol': substructure, 'atom_indices': atom_indices}


def visualize_mcs_and_non_common_parts(mol1_file, mol2_file, output_dir=".", only_non_common_rings=False, separate_rings=False):
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
    only_non_common_rings : bool, optional
        Trueの場合、複合環のうち非共通な環のみをプローブ化する, by default False
    separate_rings : bool, optional
        Trueの場合、非共通な環を別々のプローブとして抽出する, by default False
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
    non_common_parts = extract_non_common_parts(mol1_file, mol2_file, only_non_common_rings)
    
    # 結果を保存
    os.makedirs(output_dir, exist_ok=True)
    
    # 非共通部分を保存
    if separate_rings and only_non_common_rings:
        # 複数のサブ構造を保存
        for i, substructure in enumerate(non_common_parts['mol1']):
            writer = Chem.SDWriter(os.path.join(output_dir, f"non_common_mol1_ring{i+1}.sdf"))
            writer.write(substructure['mol'])
            writer.close()
        
        for i, substructure in enumerate(non_common_parts['mol2']):
            writer = Chem.SDWriter(os.path.join(output_dir, f"non_common_mol2_ring{i+1}.sdf"))
            writer.write(substructure['mol'])
            writer.close()
        
        print(f"非共通部分を {output_dir} に保存しました")
        print(f"  mol1の非共通部分: non_common_mol1_ring*.sdf")
        print(f"  mol2の非共通部分: non_common_mol2_ring*.sdf")
    else:
        # 単一のサブ構造を保存
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
    parser.add_argument('--only-non-common-rings', action='store_true',
                        help='複合環のうち非共通な環のみをプローブ化する')
    parser.add_argument('--separate-rings', action='store_true',
                        help='非共通な環を別々のプローブとして抽出する')
    args = parser.parse_args()
    
    # MCSと非共通部分を可視化
    visualize_mcs_and_non_common_parts(args.mol1, args.mol2, args.output_dir, args.only_non_common_rings, args.separate_rings)
    
    # 非共通部分を抽出
    non_common_parts = extract_non_common_parts(args.mol1, args.mol2, args.only_non_common_rings, args.separate_rings)
    
    # 非共通部分の原子数を表示
    if args.separate_rings and args.only_non_common_rings:
        # 複数のサブ構造の場合
        print("mol1の非共通部分の環ごとの原子数:")
        for i, substructure in enumerate(non_common_parts['mol1']):
            print(f"  環{i+1}: {len(substructure['atom_indices'])}原子")
        
        print("mol2の非共通部分の環ごとの原子数:")
        for i, substructure in enumerate(non_common_parts['mol2']):
            print(f"  環{i+1}: {len(substructure['atom_indices'])}原子")
    else:
        # 単一のサブ構造の場合
        print(f"mol1の非共通部分の原子数: {len(non_common_parts['mol1']['atom_indices'])}")
        print(f"mol2の非共通部分の原子数: {len(non_common_parts['mol2']['atom_indices'])}")


if __name__ == "__main__":
    main()