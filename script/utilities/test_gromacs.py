import pytest
import numpy as np
import warnings
from pathlib import Path

from script.utilities.gromacs import GroAtom, Gro, ATOMIC_NUMBER, ATOMIC_WEIGHT
from scipy import constants

class TestGroAtom:
    @pytest.fixture
    def empty_atom(self):
        """空のGroAtomインスタンスを提供するフィクスチャ"""
        return GroAtom()

    @pytest.fixture
    def normal_atom(self):
        """通常の原子情報を持つGroAtomインスタンスを提供するフィクスチャ"""
        return GroAtom("    1WAT     OW    1   0.000   0.000   0.000")

    def test_empty_initialization(self, empty_atom):
        """空のGroAtomインスタンスの作成テスト
        
        期待される動作:
        - 全ての属性が適切なデフォルト値で初期化されること
        """
        assert empty_atom.resi == -1
        assert empty_atom.resn == ""
        assert empty_atom.atomtype == ""
        assert empty_atom.atom_id == -1
        assert empty_atom.point == pytest.approx(np.zeros(3))
        assert empty_atom.velocity == pytest.approx(np.zeros(3))
        assert empty_atom.comment == ""
        assert empty_atom.atomic_mass == pytest.approx(0.0)

    @pytest.mark.parametrize("atom_str, expected", [
        (
            "    1WAT     OW    1   0.000   0.000   0.000",
            {
                "resi": 1,
                "resn": "WAT",
                "atomtype": "OW",
                "atom_id": 1,
                "point": np.array([0.0, 0.0, 0.0]),
                "velocity": np.array([0.0, 0.0, 0.0]),
                "comment": "",
                "atomic_mass": ATOMIC_WEIGHT[ATOMIC_NUMBER["O"]]
            }
        ),
        (
            "    1WAT     OW    1   0.000   0.000   0.000 ; water oxygen",
            {
                "resi": 1,
                "resn": "WAT",
                "atomtype": "OW",
                "atom_id": 1,
                "point": np.array([0.0, 0.0, 0.0]),
                "velocity": np.array([0.0, 0.0, 0.0]),
                "comment": "water oxygen",
                "atomic_mass": ATOMIC_WEIGHT[ATOMIC_NUMBER["O"]]
            }
        ),
        (
            "    1WAT     OW    1   0.000   0.000   0.000   1.000   2.000   3.000",
            {
                "resi": 1,
                "resn": "WAT",
                "atomtype": "OW",
                "atom_id": 1,
                "point": np.array([0.0, 0.0, 0.0]),
                "velocity": np.array([1.0, 2.0, 3.0]),
                "comment": "",
                "atomic_mass": ATOMIC_WEIGHT[ATOMIC_NUMBER["O"]]
            }
        )
    ])
    def test_parse_atom_variations(self, atom_str, expected):
        """様々な形式の原子情報の解析テスト
        
        パラメータ:
        - atom_str: 解析する原子情報文字列
        - expected: 期待される属性値の辞書
        
        期待される動作:
        - 全ての属性が期待される値で正しく設定されること
        """
        atom = GroAtom(atom_str)
        assert atom.resi == expected["resi"]
        assert atom.resn == expected["resn"]
        assert atom.atomtype == expected["atomtype"]
        assert atom.atom_id == expected["atom_id"]
        np.testing.assert_array_almost_equal(atom.point, expected["point"])
        np.testing.assert_array_almost_equal(atom.velocity, expected["velocity"])
        assert atom.comment == expected["comment"]
        assert atom.atomic_mass == pytest.approx(expected["atomic_mass"])

    @pytest.mark.parametrize("invalid_str", [
        "    1WAT     OW    1   0.000   0.000",  # 座標不足
        "    1WAT     OW    1",  # 座標なし
    ])
    def test_parse_invalid_coordinates(self, invalid_str):
        """無効な座標情報の解析テスト
        
        パラメータ:
        - invalid_str: 無効な原子情報文字列
        
        期待される動作:
        - RuntimeErrorが発生すること
        """
        with pytest.raises(RuntimeError, match="the dimension of atom coordinates/velocities are wrong"):
            GroAtom(invalid_str)

    def test_unknown_atomtype(self):
        """未知の原子タイプの処理テスト
        
        期待される動作:
        - RuntimeWarningが発生すること
        - atomic_massがデフォルト値に設定されること
        """
        warning_message = "atomtype XX is not matched to any atom names"
        with pytest.warns(RuntimeWarning, match=warning_message):
            atom = GroAtom("    1WAT     XX    1   0.000   0.000   0.000")
            assert atom.atomic_mass == pytest.approx(ATOMIC_WEIGHT[-1])

    def test_string_representation(self, normal_atom):
        """文字列表現の生成テスト
        
        期待される動作:
        - 正しいフォーマットの文字列が生成されること
        """
        expected = "    1  WAT   OW    1   0.000   0.000   0.000"
        assert str(normal_atom) == expected

class TestGro:
    @pytest.fixture
    def gro_content(self):
        """テスト用のGROファイル内容を提供するフィクスチャ"""
        return """Simple water system
3
    1WAT     OW    1   0.000   0.000   0.000
    1WAT    HW1    2   0.100   0.000   0.000
    1WAT    HW2    3   0.000   0.100   0.000
   5.0   5.0   5.0
"""

    @pytest.fixture
    def mock_gro_file(self, gro_content, tmp_path):
        """一時的なGROファイルを作成するフィクスチャ"""
        gro_file = tmp_path / "test.gro"
        gro_file.write_text(gro_content)
        return gro_file

    def test_empty_initialization(self):
        """空のGroインスタンスの作成テスト"""
        gro = Gro()
        assert gro.description == ""
        assert gro.natoms == 0
        assert gro.box_size == [0, 0, 0]
        assert len(gro.atoms) == 0

    def test_file_parsing(self, mock_gro_file):
        """GROファイルの解析テスト"""
        gro = Gro(str(mock_gro_file))
        assert gro.description == "Simple water system"
        assert gro.natoms == 3
        assert gro.box_size == [5.0, 5.0, 5.0]
        assert len(gro.atoms) == 3
        assert [atom.atomtype for atom in gro.atoms] == ["OW", "HW1", "HW2"]

    @pytest.mark.parametrize("search_params, expected_count", [
        ({"resi": 1}, 3),
        ({"resn": "WAT"}, 3),
        ({"atomtype": "OW"}, 1),
        ({"atom_id": 2}, 1),
        ({"resi": 1, "atomtype": "OW"}, 1),
    ])
    def test_get_atoms(self, mock_gro_file, search_params, expected_count):
        """原子検索機能のテスト"""
        gro = Gro(str(mock_gro_file))
        atoms = gro.get_atoms(**search_params)
        assert len(atoms) == expected_count

    def test_add_atom(self, mock_gro_file):
        """原子追加機能のテスト"""
        gro = Gro(str(mock_gro_file))
        new_atom = GroAtom("    2  WAT   OW    4   1.000   1.000   1.000")
        gro.add_atom(new_atom)
        
        assert gro.natoms == 4
        assert gro.atoms[-1].resi == 2
        assert gro.atoms[-1].atom_id == 4

    def test_add_invalid_atom(self):
        """無効な原子追加のテスト"""
        gro = Gro()
        with pytest.raises(TypeError, match="the input is NON-GRO_ATOM"):
            gro.add_atom("not a GroAtom")

    def test_molar_concentration(self, mock_gro_file):
        """モル濃度計算のテスト"""
        gro = Gro(str(mock_gro_file))
        molar = gro.molar("WAT")
        # ボックスサイズ 5.0 nm × 5.0 nm × 5.0 nm = 125 nm3
        expected_molar = (1 / constants.N_A) / ((5*constants.nano)**3 * 10**3)  # nm3 -> cm3 * 10^3 = L
        assert molar == pytest.approx(expected_molar, rel=1e-5)

    def test_string_representation(self, mock_gro_file):
        """文字列表現の生成テスト"""
        gro = Gro(str(mock_gro_file))
        expected = """Simple water system
     3
    1  WAT   OW    1   0.000   0.000   0.000
    1  WAT  HW1    2   0.100   0.000   0.000
    1  WAT  HW2    3   0.000   0.100   0.000
   5.00000     5.00000     5.00000
"""
        assert str(gro) == expected