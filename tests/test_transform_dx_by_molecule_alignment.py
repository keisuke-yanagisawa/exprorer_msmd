#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
script.analysis.inverse_msmdモジュールのtransform_dx_based_on_molecule_alignmentのテストコード
"""

import os
from pathlib import Path

import gridData
import pytest

# 関数をscript.analysis.inverse_msmdからインポート
from script.analysis.inverse_msmd import transform_dx_based_on_molecule_alignment


def test_transform_dx_by_molecule_alignment_creates_valid_output():
    """
    transform_dx_based_on_molecule_alignment関数が有効な出力を生成することをテストする
    """
    # テスト用のファイルパス
    data_dir = Path(__file__).parent / "data"
    compound1_file = data_dir / "4HW2A_lig.sdf"
    compound2_file = data_dir / "4HW3A_probe.sdf"
    dx_file = data_dir / "4HW3A_probe_ALA_profile.dx"

    # 一時的な出力ファイル
    output_dx = data_dir / "output_test.dx"

    # 分子の読み込み
    from rdkit import Chem

    mol1 = Chem.SDMolSupplier(str(compound1_file))[0]
    mol2 = Chem.SDMolSupplier(str(compound2_file))[0]

    # dxファイルの読み込み
    grid = gridData.Grid(dx_file)

    # 関数を実行
    result = transform_dx_based_on_molecule_alignment(mol1, mol2, grid, output_dx)

    # 出力ファイルが存在することを確認
    assert os.path.exists(output_dx)

    # 出力ファイルが読み込めることを確認
    grid = gridData.Grid(output_dx)

    # グリッドデータが有効であることを確認
    assert grid.grid.shape[0] > 0
    assert grid.grid.shape[1] > 0
    assert grid.grid.shape[2] > 0

    # 後片付け
    if os.path.exists(output_dx):
        os.remove(output_dx)
