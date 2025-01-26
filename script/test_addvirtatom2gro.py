import pytest
import numpy as np
from pathlib import Path

from script.utilities import gromacs
from script.addvirtatom2gro import center_of_mass, addvirtatom2gro


@pytest.fixture
def single_carbon_atom():
    """単一の炭素原子のテストデータ"""
    atom = gromacs.GroAtom()
    atom.point = np.array([1.0, 1.0, 1.0])
    atom.atomic_mass = 12.0  # 炭素原子
    return atom


@pytest.fixture
def water_molecule():
    """水分子（O + 2H）のテストデータ"""
    oxygen = gromacs.GroAtom()
    oxygen.point = np.array([0.0, 0.0, 0.0])
    oxygen.atomic_mass = 16.0  # 酸素原子

    hydrogen1 = gromacs.GroAtom()
    hydrogen1.point = np.array([0.1, 0.0, 0.0])
    hydrogen1.atomic_mass = 1.0  # 水素原子

    hydrogen2 = gromacs.GroAtom()
    hydrogen2.point = np.array([0.0, 0.1, 0.0])
    hydrogen2.atomic_mass = 1.0  # 水素原子

    return [oxygen, hydrogen1, hydrogen2]


@pytest.fixture
def co_molecule():
    """一酸化炭素分子（C + O）のテストデータ"""
    carbon = gromacs.GroAtom()
    carbon.point = np.array([0.0, 0.0, 0.0])
    carbon.atomic_mass = 12.0  # 炭素原子

    oxygen = gromacs.GroAtom()
    oxygen.point = np.array([1.0, 0.0, 0.0])
    oxygen.atomic_mass = 16.0  # 酸素原子

    return [carbon, oxygen]


def test_single_atom(single_carbon_atom):
    """単一原子の質量中心計算テスト"""
    result = center_of_mass([single_carbon_atom])
    np.testing.assert_array_almost_equal(result, single_carbon_atom.point)


def test_multiple_atoms(water_molecule):
    """複数原子（水分子）の質量中心計算テスト"""
    result = center_of_mass(water_molecule)
    # 水分子の質量中心を手動で計算
    expected = (16.0 * water_molecule[0].point +
               1.0 * water_molecule[1].point +
               1.0 * water_molecule[2].point) / 18.0
    np.testing.assert_array_almost_equal(result, expected)


def test_different_masses(co_molecule):
    """異なる質量を持つ原子の質量中心計算テスト"""
    result = center_of_mass(co_molecule)
    # 質量中心は重い方（酸素）に寄る
    expected = np.array([16.0 / (12.0 + 16.0), 0.0, 0.0])
    np.testing.assert_array_almost_equal(result, expected)


@pytest.fixture
def basic_gro():
    """単一の水分子を含むGROファイル"""
    return """Simple water system
3
    1WAT     OW    1   0.000   0.000   0.000
    1WAT    HW1    2   0.100   0.000   0.000
    1WAT    HW2    3   0.000   0.100   0.000
   5.0   5.0   5.0
"""


@pytest.fixture
def multi_probe_gro():
    """複数の水分子を含むGROファイル"""
    return """System with multiple probes
6
    1WAT     OW    1   0.000   0.000   0.000
    1WAT    HW1    2   0.100   0.000   0.000
    1WAT    HW2    3   0.000   0.100   0.000
    2WAT     OW    4   1.000   1.000   1.000
    2WAT    HW1    5   1.100   1.000   1.000
    2WAT    HW2    6   1.000   1.100   1.000
   5.0   5.0   5.0
"""


@pytest.fixture
def empty_gro():
    """空のGROファイル"""
    return """Empty system
0
5.00000     5.00000     5.00000
"""


def test_basic_conversion(basic_gro):
    """基本的なGROファイルの変換テスト"""
    result = addvirtatom2gro(basic_gro, "WAT")
    
    # 仮想原子が追加されていることを確認
    assert "VIS" in result
    # 仮想原子の座標が質量中心に近いことを確認
    assert "0.006   0.006   0.000" in result  # 水分子の質量中心


def test_multiple_probes(multi_probe_gro):
    """複数のプローブ分子を含むGROファイルの変換テスト"""
    result = addvirtatom2gro(multi_probe_gro, "WAT")
    
    # 2つの仮想原子が追加されていることを確認
    vis_count = result.count("VIS")
    assert vis_count == 2
    
    # 両方の水分子の質量中心に仮想原子が配置されていることを確認
    assert "0.006   0.006   0.000" in result  # 1番目の水分子
    assert "1.006   1.006   1.000" in result  # 2番目の水分子


def test_no_probe(basic_gro):
    """プローブが存在しない場合のテスト"""
    result = addvirtatom2gro(basic_gro, "ETH")  # 存在しないプローブID
    
    # 仮想原子が追加されていないことを確認
    assert "VIS" not in result
    # 元のGROファイルの内容が保持されていることを確認
    assert len(result.split("\n")) == len(basic_gro.split("\n"))


def test_empty_gro(empty_gro):
    """空のGROファイルの処理テスト"""
    result = addvirtatom2gro(empty_gro, "WAT")
    
    # 期待値と実際の出力を正規化して比較
    def normalize_gro(gro_str):
        lines = gro_str.strip().split("\n")
        return "\n".join([
            lines[0],  # システム名
            lines[1].strip(),  # 原子数
            lines[2].strip()  # ボックスサイズ
        ])
    assert normalize_gro(result) == normalize_gro(empty_gro)