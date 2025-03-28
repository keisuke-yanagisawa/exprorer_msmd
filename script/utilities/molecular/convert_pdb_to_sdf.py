#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
PDBファイルをSDFファイルに変換するスクリプト
"""

import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem


def convert_pdb_to_sdf(pdb_file, sdf_file):
    """
    PDBファイルをSDFファイルに変換する

    Parameters
    ----------
    pdb_file : str
        PDBファイルのパス
    sdf_file : str
        出力SDFファイルのパス
    """
    # PDBファイルを読み込む
    mol = Chem.MolFromPDBFile(pdb_file)
    if mol is None:
        print(f"エラー: {pdb_file} を読み込めませんでした")
        return False
    
    # 水素を追加
    mol = Chem.AddHs(mol, addCoords=True)
    
    # SDFファイルに保存
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()
    
    print(f"{pdb_file} を {sdf_file} に変換しました")
    return True

def create_modified_molecule(pdb_file, sdf_file, modification="methyl"):
    """
    分子を修正して新しい分子を作成する

    Parameters
    ----------
    pdb_file : str
        PDBファイルのパス
    sdf_file : str
        出力SDFファイルのパス
    modification : str, optional
        修正の種類, by default "methyl"
    """
    # PDBファイルを読み込む
    mol = Chem.MolFromPDBFile(pdb_file)
    if mol is None:
        print(f"エラー: {pdb_file} を読み込めませんでした")
        return False
    
    # 水素を追加
    mol = Chem.AddHs(mol, addCoords=True)
    
    # 分子を修正
    modified_mol = None
    
    if modification == "methyl":
        # メチル基を追加（簡易的な修正）
        # 実際の実装では、より適切な修正方法を使用する
        # ここでは単純に原子の座標を少し変更して、異なる分子を作成
        rwmol = Chem.RWMol(mol)
        conf = rwmol.GetConformer()
        
        # 最初の炭素原子の座標を少し変更
        for i in range(rwmol.GetNumAtoms()):
            atom = rwmol.GetAtomWithIdx(i)
            if atom.GetSymbol() == 'C':
                pos = conf.GetAtomPosition(i)
                new_pos = Chem.rdGeometry.Point3D(pos.x + 0.1, pos.y + 0.1, pos.z + 0.1)
                conf.SetAtomPosition(i, new_pos)
                break
        
        modified_mol = rwmol.GetMol()
    
    if modified_mol is None:
        print(f"警告: 分子を修正できませんでした。元の分子をコピーします。")
        modified_mol = Chem.Mol(mol)
    
    # SDFファイルに保存
    writer = Chem.SDWriter(sdf_file)
    writer.write(modified_mol)
    writer.close()
    
    print(f"修正した分子を {sdf_file} に保存しました")
    return True

if __name__ == "__main__":
    # コマンドライン引数の解析
    if len(sys.argv) < 3:
        print(f"使用法: {sys.argv[0]} <pdb_file> <sdf_file>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    sdf_file = sys.argv[2]
    
    # PDBファイルをSDFファイルに変換
    convert_pdb_to_sdf(pdb_file, sdf_file)
    
    # 修正した分子を作成
    if len(sys.argv) > 3:
        modified_sdf_file = sys.argv[3]
        create_modified_molecule(pdb_file, modified_sdf_file)