#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
SMILES表記からSDFファイルを生成するスクリプト
"""

import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_sdf(smiles, output_file, name=None):
    """
    SMILES表記からSDFファイルを生成する

    Parameters
    ----------
    smiles : str
        SMILES表記
    output_file : str
        出力SDFファイルパス
    name : str, optional
        分子名, by default None
    """
    # SMILES表記から分子を生成
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"エラー: SMILES表記 {smiles} から分子を生成できませんでした")
        return False
    
    # 分子名を設定
    if name:
        mol.SetProp("_Name", name)
    
    # 水素を追加
    mol = Chem.AddHs(mol)
    
    # 3D座標を生成
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    
    # SDFファイルに保存
    writer = Chem.SDWriter(output_file)
    writer.write(mol)
    writer.close()
    
    print(f"SMILES表記 {smiles} からSDFファイル {output_file} を生成しました")
    return True

def main():
    """
    メイン関数
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='SMILES表記からSDFファイルを生成するスクリプト')
    parser.add_argument('--smiles', required=True, help='SMILES表記')
    parser.add_argument('--output', required=True, help='出力SDFファイルパス')
    parser.add_argument('--name', help='分子名')
    args = parser.parse_args()
    
    smiles_to_sdf(args.smiles, args.output, args.name)

if __name__ == "__main__":
    main()