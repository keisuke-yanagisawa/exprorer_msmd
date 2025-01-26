import os
import pytest
from pathlib import Path
from typing import Literal, cast

from script.generate_msmd_system import (
    _create_frcmod,
    calculate_boxsize,
    generate_msmd_system,
    protein_pdb_preparation
)
from script.setting import parse_yaml

# テストデータのパスを設定
TEST_DATA_DIR = Path("script/test_data")

@pytest.fixture
def test_files():
    return {
        'rst7': TEST_DATA_DIR / "tripeptide.rst7",
        'pdb': TEST_DATA_DIR / "tripeptide.pdb",
        'setting': TEST_DATA_DIR / "setting.yaml",
        'expected_parm7': TEST_DATA_DIR / "tripeptide_A11.parm7",
        'expected_rst7': TEST_DATA_DIR / "tripeptide_A11.rst7",
        'mol2': TEST_DATA_DIR / "A11.mol2",
        'expected_frcmod': TEST_DATA_DIR / "A11.frcmod"
    }

def compare_file_contents(file1: Path, file2: Path, skip_first_line: bool = False) -> None:
    """2つのファイルの内容を比較する

    Args:
        file1: 比較対象ファイル1
        file2: 比較対象ファイル2
        skip_first_line: 最初の行をスキップするかどうか
    """
    try:
        with open(file1) as f1, open(file2) as f2:
            content1 = f1.readlines()[1:] if skip_first_line else f1.readlines()
            content2 = f2.readlines()[1:] if skip_first_line else f2.readlines()
            
            assert content1 == content2, (
                f"File contents do not match:\n"
                f"File 1 ({file1}): {len(content1)} lines\n"
                f"File 2 ({file2}): {len(content2)} lines"
            )
    except FileNotFoundError as e:
        raise AssertionError(f"File not found: {e.filename}")
    except Exception as e:
        raise AssertionError(f"Error comparing files: {str(e)}")

class TestBoxSizeCalculation:
    def test_calculate_boxsize(self, test_files):
        box_size = calculate_boxsize(test_files['rst7'])
        expected = 16.4927710
        assert box_size == pytest.approx(expected, rel=1e-6)

    def test_calculate_boxsize_error(self, test_files):
        """PDBファイルはrst7形式ではないため、ボックスサイズの計算に失敗することを確認"""
        with pytest.raises(ValueError):
            calculate_boxsize(test_files['pdb'])

class TestGenerateMsmdSystem:
    def test_generate_msmd_system_parm7(self, test_files):
        """生成されたparm7ファイルが期待通りであることを確認"""
        settings = parse_yaml(test_files['setting'])
        parm7, _ = generate_msmd_system(settings, seed=1)
        compare_file_contents(parm7, test_files['expected_parm7'], skip_first_line=True)

    def test_generate_msmd_system_rst7(self, test_files):
        """生成されたrst7ファイルが期待通りであることを確認"""
        settings = parse_yaml(test_files['setting'])
        _, rst7 = generate_msmd_system(settings, seed=1)
        compare_file_contents(rst7, test_files['expected_rst7'])

    def test_protein_pdb_preparation(self, test_files, tmp_path):
        # テスト用の一時的なPDBファイルを作成
        test_pdb = tmp_path / "test.pdb"
        with open(test_pdb, "w") as f:
            f.write("ATOM      1  N   ALA A   1      27.409  24.354   9.020  1.00  0.00           N  \n")
            f.write("ATOM      2  OXT ALA A   1      27.409  24.354   9.020  1.00  0.00           O  \n")
            f.write("ANISOU    1  N   ALA A   1     2406   2304   2278    -28     54   -167       N  \n")

        result_path = protein_pdb_preparation(test_pdb)
        
        # 結果の検証
        with open(result_path) as f:
            lines = f.readlines()
        
        assert len(lines) == 1
        assert "OXT" not in lines[0]
        assert "ANISOU" not in lines[0]

class TestCreateFrcmod:
    def test_create_frcmod(self, test_files):
        """frcmodファイルが正しく生成されることを確認"""
        atomtype: Literal["gaff", "gaff2"] = "gaff2"
        frcmod_path = _create_frcmod(test_files['mol2'], atomtype)
        compare_file_contents(frcmod_path, test_files['expected_frcmod'])

    def test_invalid_atomtype(self, test_files):
        with pytest.raises(ValueError):
            invalid_type = cast(Literal["gaff", "gaff2"], "invalid")
            _create_frcmod(test_files['mol2'], invalid_type)

    def test_mol2file_does_not_exist(self):
        atomtype: Literal["gaff", "gaff2"] = "gaff2"
        with pytest.raises(FileNotFoundError):
            _create_frcmod(Path("INVALID"), atomtype)

    def test_invalid_mol2file(self, test_files):
        atomtype: Literal["gaff", "gaff2"] = "gaff2"
        with pytest.raises(ValueError):
            _create_frcmod(test_files['pdb'], atomtype)