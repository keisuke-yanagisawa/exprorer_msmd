#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
分子の特定の原子を別の原子に置換するシンプルなスクリプト
"""

import os

from rdkit import Chem
from rdkit.Chem import AllChem


def replace_atom_in_molecule(mol_file, atom_idx, new_atom_symbol, output_dir="output"):
    """
    分子の特定の原子を別の原子に置換する
    
    Parameters
    ----------
    mol_file : str
        分子のSDFファイルパス
    atom_idx : int
        置換対象原子のインデックス
    new_atom_symbol : str
        新しい原子のシンボル（例: 'N', 'O'）
    output_dir : str, optional
        出力ディレクトリ, by default "output"
    
    Returns
    -------
    tuple
        (元の分子のファイルパス, 置換後の分子のファイルパス)
    """
    # 出力ディレクトリを作成
    os.makedirs(output_dir, exist_ok=True)
    
    # 分子の読み込み
    mol_supplier = Chem.SDMolSupplier(mol_file)
    mol = mol_supplier[0]
    
    if mol is None:
        raise ValueError(f"分子ファイル {mol_file} を読み込めませんでした")
    
    # 元の分子を保存
    original_file = os.path.join(output_dir, "original_molecule.sdf")
    writer = Chem.SDWriter(original_file)
    writer.write(mol)
    writer.close()
    
    print(f"元の分子を {original_file} に保存しました")
    
    # 原子を置換
    # 新しい分子を作成
    rwmol = Chem.RWMol(mol)
    
    # 置換対象の原子を取得
    atom = rwmol.GetAtomWithIdx(atom_idx)
    
    # 原子の元素を変更
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
    
    # 原子の元素を変更
    atom.SetAtomicNum(atomic_num)
    
    # 分子を更新
    replaced_mol = rwmol.GetMol()
    
    # 分子を正規化
    try:
        Chem.SanitizeMol(replaced_mol)
    except Exception as e:
        print(f"警告: 分子の正規化に失敗しました: {e}")
    
    # 置換後の分子を保存
    replaced_file = os.path.join(output_dir, "replaced_molecule.sdf")
    writer = Chem.SDWriter(replaced_file)
    writer.write(replaced_mol)
    writer.close()
    
    print(f"置換後の分子を {replaced_file} に保存しました")
    
    return original_file, replaced_file


def main():
    """
    メイン関数
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='分子の特定の原子を別の原子に置換する')
    parser.add_argument('--mol', required=True, help='分子のSDFファイルパス')
    parser.add_argument('--atom-idx', type=int, required=True, help='置換対象原子のインデックス')
    parser.add_argument('--new-atom', required=True, help='新しい原子のシンボル（例: N, O）')
    parser.add_argument('--output-dir', default="output", help='出力ディレクトリ')
    args = parser.parse_args()
    
    # 分子の処理
    original_file, replaced_file = replace_atom_in_molecule(
        args.mol, args.atom_idx, args.new_atom, args.output_dir
    )
    
    print("処理が完了しました")
    print(f"元の分子: {original_file}")
    print(f"置換後の分子: {replaced_file}")


if __name__ == "__main__":
    main()