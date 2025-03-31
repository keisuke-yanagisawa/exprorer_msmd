#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
script.analysis.inverse_msmdモジュールのextract_non_common_partsのテストコード
"""

from pathlib import Path

import pytest
from rdkit import Chem

# 関数をscript.analysis.inverse_msmdからインポート
from script.analysis.inverse_msmd import extract_non_common_parts


def test_extract_non_common_parts_returns_valid_result():
    """
    extract_non_common_parts関数が有効な結果を返すことをテストする
    """
    # テスト用のファイルパス
    data_dir = Path(__file__).parent / "data"
    mol1_file = data_dir / "4HW2A_lig.sdf"
    mol2_file = data_dir / "4HW3A_lig.sdf"

    # 分子の読み込み
    mol1 = Chem.SDMolSupplier(str(mol1_file))[0]
    mol2 = Chem.SDMolSupplier(str(mol2_file))[0]

    # 関数を実行
    result = extract_non_common_parts(mol1, mol2)

    # 結果が辞書であることを確認
    assert isinstance(result, dict)
    # 結果に "mol1" と "mol2" のキーが含まれていることを確認
    assert "mol1" in result
    assert "mol2" in result
    # 結果の各要素がリストであることを確認
    assert isinstance(result["mol1"], list)
    assert isinstance(result["mol2"], list)
