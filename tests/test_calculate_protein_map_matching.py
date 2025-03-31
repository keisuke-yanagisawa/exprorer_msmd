#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
script.analysis.inverse_msmdモジュールのcalculate_protein_map_matchingのテストコード
"""

from pathlib import Path

import gridData
import pytest

# 関数をscript.analysis.inverse_msmdからインポート
from script.analysis.inverse_msmd import calculate_protein_map_matching
from script.utilities.Bio.PDB import get_structure


def test_calculate_protein_map_matching_returns_valid_score():
    """
    calculate_protein_map_matching関数が有効なスコアを返すことをテストする
    """
    # テスト用のファイルパス
    data_dir = Path(__file__).parent / "data"
    protein_file = data_dir / "4HW2A_protein.pdb"
    profile_file = data_dir / "4HW3A_probe_ALA_profile.dx"

    # プロファイルファイルの辞書を作成
    profile_files = {"ALA": gridData.Grid(profile_file)}

    # 構造の取得
    protein = get_structure(protein_file)

    # 関数を実行
    score = calculate_protein_map_matching(protein, profile_files)

    # スコアが数値であることを確認
    assert isinstance(score, (int, float))
    # スコアが有限であることを確認
    assert not (score == float("inf") or score == float("-inf") or score != score)
