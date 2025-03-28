#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ZINC000000330081.sdfの分子について、4番目のCをNに置き換えた場合の
結合親和性の変化をquantitative inverse MSMDで見積もるためのスクリプト
"""

import os
import sys

from substructure_extractor import process_molecule_for_msmd


def main():
    """
    メイン関数
    """
    # 入力ファイル
    mol_file = "example/ZINC000000330081.sdf"
    
    # 置換対象原子のインデックス（4番目のC）
    target_atom_idx = 3  # 0-indexedなので、4番目は3
    
    # 新しい原子のシンボル
    new_atom_symbol = "N"
    
    print(f"処理を開始します: {mol_file}")
    print(f"置換対象原子: インデックス {target_atom_idx} (4番目のC)")
    print(f"置換後の原子: {new_atom_symbol}")
    
    try:
        # 分子の処理
        original_file, replaced_file = process_molecule_for_msmd(
            mol_file, target_atom_idx, new_atom_symbol
        )
        
        print("\n処理が完了しました")
        print(f"元の部分構造: {original_file}")
        print(f"置換後の部分構造: {replaced_file}")
        print("\n次のステップ: 生成された部分構造に対してquantitative inverse MSMDを実行してください")
        
    except Exception as e:
        print(f"エラーが発生しました: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()