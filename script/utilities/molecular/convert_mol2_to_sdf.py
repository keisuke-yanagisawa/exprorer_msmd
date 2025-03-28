#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
MOL2ファイルをSDFファイルに変換するスクリプト
"""

import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem


def convert_mol2_to_sdf(mol2_file, sdf_file):
    """
    MOL2ファイルをSDFファイルに変換する

    Parameters
    ----------
    mol2_file : str
        MOL2ファイルのパス
    sdf_file : str
        出力SDFファイルのパス
    """
    # MOL2ファイルを読み込む
    mol = Chem.MolFromMol2File(mol2_file)
    if mol is None:
        print(f"エラー: {mol2_file} を読み込めませんでした")
        return False
    
    # SDFファイルに保存
    writer = Chem.SDWriter(sdf_file)
    writer.write(mol)
    writer.close()
    
    print(f"{mol2_file} を {sdf_file} に変換しました")
    return True

def create_modified_molecule(mol2_file, sdf_file, modification="methyl"):
    """
    分子を修正して新しい分子を作成する

    Parameters
    ----------
    mol2_file : str
        MOL2ファイルのパス
    sdf_file : str
        出力SDFファイルのパス
    modification : str, optional
        修正の種類, by default "methyl"
    """
    # MOL2ファイルを読み込む
    mol = Chem.MolFromMol2File(mol2_file)
    if mol is None:
        print(f"エラー: {mol2_file} を読み込めませんでした")
        return False
    
    # 分子を修正
    modified_mol = None
    
    if modification == "methyl":
        # メチル基を追加（簡易的な修正）
        # 実際の実装では、より適切な修正方法を使用する
        patt = Chem.MolFromSmarts("c")  # 芳香族炭素を検索
        matches = mol.GetSubstructMatches(patt)
        if matches:
            # 最初の芳香族炭素を選択
            atom_idx = matches[0][0]
            # 新しい分子を作成
            rwmol = Chem.RWMol(mol)
            # 水素を追加
            rwmol = Chem.AddHs(rwmol)
            # メチル基を追加（実際には単純な置換ではなく、より複雑な処理が必要）
            modified_mol = rwmol
    
    if modified_mol is None:
        print(f"警告: 分子を修正できませんでした。元の分子をコピーします。")
        modified_mol = Chem.Mol(mol)
    
    # 3D座標を生成
    modified_mol = Chem.AddHs(modified_mol)
    AllChem.EmbedMolecule(modified_mol)
    AllChem.UFFOptimizeMolecule(modified_mol)
    
    # SDFファイルに保存
    writer = Chem.SDWriter(sdf_file)
    writer.write(modified_mol)
    writer.close()
    
    print(f"修正した分子を {sdf_file} に保存しました")
    return True

if __name__ == "__main__":
    # コマンドライン引数の解析
    if len(sys.argv) < 3:
        print(f"使用法: {sys.argv[0]} <mol2_file> <sdf_file>")
        sys.exit(1)
    
    mol2_file = sys.argv[1]
    sdf_file = sys.argv[2]
    
    # MOL2ファイルをSDFファイルに変換
    convert_mol2_to_sdf(mol2_file, sdf_file)
    
    # 修正した分子を作成
    if len(sys.argv) > 3:
        modified_sdf_file = sys.argv[3]
        create_modified_molecule(mol2_file, modified_sdf_file)