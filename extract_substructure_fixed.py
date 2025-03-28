#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
分子から部分構造を抽出するスクリプト（修正版）
"""

import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem


def extract_substructure_around_atom(mol_file, atom_idx, radius=2, output_dir="output"):
    """
    指定した原子の周りの部分構造を抽出する
    
    Parameters
    ----------
    mol_file : str
        分子のSDFファイルパス
    atom_idx : int
        中心となる原子のインデックス
    radius : int, optional
        取り込む周辺原子の範囲（結合距離）, by default 2
    output_dir : str, optional
        出力ディレクトリ, by default "output"
    
    Returns
    -------
    str
        抽出された部分構造のSDFファイルパス
    """
    # 出力ディレクトリを作成
    os.makedirs(output_dir, exist_ok=True)
    
    # 分子の読み込み
    mol_supplier = Chem.SDMolSupplier(mol_file)
    mol = mol_supplier[0]
    
    if mol is None:
        raise ValueError(f"分子ファイル {mol_file} を読み込めませんでした")
    
    # 中心原子から指定した半径内の原子を取得
    atoms = set()
    atoms.add(atom_idx)
    
    # 指定した半径まで原子を追加
    current_atoms = {atom_idx}
    for r in range(radius):
        next_atoms = set()
        for atom_id in current_atoms:
            atom = mol.GetAtomWithIdx(atom_id)
            for neighbor in atom.GetNeighbors():
                next_atoms.add(neighbor.GetIdx())
        # 新しい原子を追加
        atoms.update(next_atoms)
        current_atoms = next_atoms
    
    # 環構造を保持するための処理
    # 中心原子が環構造に含まれる場合、その環全体を含める
    ring_info = mol.GetRingInfo()
    for ring_atoms in ring_info.AtomRings():
        if atom_idx in ring_atoms:
            atoms.update(ring_atoms)
    
    # 部分構造を抽出
    atoms_list = list(atoms)
    
    # 部分構造を抽出するためにRWMolを使用
    rwmol = Chem.RWMol()
    
    # 原子マッピングを作成
    atom_mapping = {}
    
    # 選択した原子を追加
    for old_idx in atoms_list:
        atom = mol.GetAtomWithIdx(old_idx)
        # 原子番号を取得
        atomic_num = atom.GetAtomicNum()
        # 新しい原子を追加
        new_idx = rwmol.AddAtom(Chem.Atom(atomic_num))
        atom_mapping[old_idx] = new_idx
        
        # 原子の特性をコピー
        new_atom = rwmol.GetAtomWithIdx(new_idx)
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetIsAromatic(atom.GetIsAromatic())
        new_atom.SetChiralTag(atom.GetChiralTag())
    
    # 結合を追加
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        if begin_idx in atoms_list and end_idx in atoms_list:
            rwmol.AddBond(
                atom_mapping[begin_idx],
                atom_mapping[end_idx],
                bond.GetBondType()
            )
    
    # RWMolから通常の分子に変換
    substructure = rwmol.GetMol()
    
    # 3D座標を設定
    if mol.GetNumConformers() > 0:
        conf = Chem.Conformer(substructure.GetNumAtoms())
        for old_idx, new_idx in atom_mapping.items():
            old_pos = mol.GetConformer().GetAtomPosition(old_idx)
            conf.SetAtomPosition(new_idx, old_pos)
        
        substructure.AddConformer(conf)
    
    # 分子を正規化（sanitize）
    try:
        Chem.SanitizeMol(substructure)
    except Exception as e:
        print(f"警告: 分子の正規化に失敗しました: {e}")
    
    # 部分構造を保存
    output_file = os.path.join(output_dir, os.path.basename(mol_file).replace(".sdf", "_substructure.sdf"))
    writer = Chem.SDWriter(output_file)
    writer.write(substructure)
    writer.close()
    
    print(f"部分構造を {output_file} に保存しました")
    
    return output_file


def main():
    """
    メイン関数
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='分子から部分構造を抽出する')
    parser.add_argument('--mol', required=True, help='分子のSDFファイルパス')
    parser.add_argument('--atom-idx', type=int, required=True, help='中心となる原子のインデックス')
    parser.add_argument('--radius', type=int, default=2, help='取り込む周辺原子の範囲（結合距離）')
    parser.add_argument('--output-dir', default="output", help='出力ディレクトリ')
    args = parser.parse_args()
    
    # 部分構造の抽出
    output_file = extract_substructure_around_atom(
        args.mol, args.atom_idx, args.radius, args.output_dir
    )
    
    print("処理が完了しました")
    print(f"部分構造: {output_file}")


if __name__ == "__main__":
    main()