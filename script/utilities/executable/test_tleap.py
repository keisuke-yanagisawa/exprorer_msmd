import os
import pytest
from pathlib import Path
from typing import Optional, List
from script.utilities.executable.tleap import TLeap
import subprocess

# テストデータの定義
A11_MOL2_CONTENT = """@<TRIPOS>MOLECULE
A11
   12    11     1     0     0
SMALL
resp


@<TRIPOS>ATOM
      1 C1           1.2730    -0.6670     0.0980 c3         1 A11    -0.551360
      2 C2           0.0000     0.0400    -0.3710 c3         1 A11     0.625940
      3 C3          -1.2730    -0.6670     0.0980 c3         1 A11    -0.551360
      4 O1           0.0000     1.4150     0.0230 oh         1 A11    -0.734700
      5 H1           1.3250    -1.6940    -0.2830 hc         1 A11     0.129000
      6 H2           1.3070    -0.7190     1.1950 hc         1 A11     0.129000
      7 H3           2.1590    -0.1230    -0.2450 hc         1 A11     0.129000
      8 H4           0.0000     0.0800    -1.4670 h1         1 A11     0.015380
      9 H5          -1.3070    -0.7190     1.1950 hc         1 A11     0.129000
     10 H6          -1.3250    -1.6940    -0.2830 hc         1 A11     0.129000
     11 H7          -2.1590    -0.1230    -0.2450 hc         1 A11     0.129000
     12 H8           0.0000     1.4340     0.9940 ho         1 A11     0.422100
@<TRIPOS>BOND
     1     1     2 1   
     2     1     5 1   
     3     1     6 1   
     4     1     7 1   
     5     2     3 1   
     6     2     4 1   
     7     2     8 1   
     8     3     9 1   
     9     3    10 1   
    10     3    11 1   
    11     4    12 1   
@<TRIPOS>SUBSTRUCTURE
     1 A11         1 TEMP              0 ****  ****    0 ROOT"""

A11_FRCMOD_CONTENT = """Remark line goes here
MASS

BOND

ANGLE

DIHE

IMPROPER

NONBON


"""

TRIPEPTIDE_PDB_CONTENT = """TITLE     2 Tri-peptide
MODEL        1
ATOM      1  N   MET     1     -23.351 -17.664   9.173  1.00 33.81           N
ATOM      2  CA  MET     1     -22.530 -17.005   8.116  1.00 32.98           C
ATOM      3  C   MET     1     -21.431 -16.158   8.740  1.00 30.36           C
ATOM      4  O   MET     1     -21.695 -15.267   9.549  1.00 29.19           O
ATOM      5  CB  MET     1     -23.410 -16.118   7.228  1.00 37.05           C
ATOM      6  CG  MET     1     -22.648 -15.341   6.150  1.00 42.54           C
ATOM      7  SD  MET     1     -23.748 -14.497   4.958  1.00 50.60           S
ATOM      8  CE  MET     1     -22.560 -13.892   3.738  1.00 48.71           C
ATOM      9  N   ASP     2     -20.194 -16.450   8.362  1.00 27.76           N
ATOM     10  CA  ASP     2     -19.045 -15.709   8.864  1.00 26.73           C
ATOM     11  C   ASP     2     -18.610 -14.660   7.837  1.00 25.55           C
ATOM     12  O   ASP     2     -19.188 -14.554   6.756  1.00 24.66           O
ATOM     13  CB  ASP     2     -17.899 -16.685   9.193  1.00 23.81           C
ATOM     14  CG  ASP     2     -17.561 -17.605   8.030  1.00 26.35           C
ATOM     15  OD1 ASP     2     -16.822 -18.593   8.225  1.00 25.32           O
ATOM     16  OD2 ASP     2     -18.036 -17.337   6.914  1.00 25.39           O
ATOM     17  N   LYS     3     -17.607 -13.868   8.190  1.00 26.45           N
ATOM     18  CA  LYS     3     -17.116 -12.826   7.295  1.00 25.92           C
ATOM     19  C   LYS     3     -16.555 -13.433   6.012  1.00 25.57           C
ATOM     20  O   LYS     3     -15.639 -14.246   6.057  1.00 25.30           O
ATOM     21  CB  LYS     3     -16.027 -12.024   8.000  1.00 27.03           C
ATOM     22  CG  LYS     3     -15.327 -10.981   7.149  1.00 25.32           C
ATOM     23  CD  LYS     3     -14.276 -10.277   7.997  1.00 29.33           C
ATOM     24  CE  LYS     3     -13.400  -9.343   7.182  1.00 30.07           C
ATOM     25  NZ  LYS     3     -14.233  -8.349   6.468  1.00 31.27           N
TER
ENDMDL"""

@pytest.fixture
def test_files(tmp_path: Path) -> dict:
    """テストに必要なファイルを一時ディレクトリに作成するフィクスチャ"""
    probe_path = tmp_path / "A11.mol2"
    frcmod_path = tmp_path / "A11.frcmod"
    box_path = tmp_path / "tripeptide.pdb"

    probe_path.write_text(A11_MOL2_CONTENT)
    frcmod_path.write_text(A11_FRCMOD_CONTENT)
    box_path.write_text(TRIPEPTIDE_PDB_CONTENT)

    return {
        "probe_path": probe_path,
        "frcmod_path": frcmod_path,
        "box_path": box_path
    }

@pytest.fixture
def tleap() -> TLeap:
    """TLeapインスタンスを提供するフィクスチャ"""
    return TLeap(debug=True)

@pytest.fixture
def tleap_with_params(tleap: TLeap, test_files: dict) -> TLeap:
    """基本パラメータが設定されたTLeapインスタンスを提供するフィクスチャ"""
    return tleap.set(
        cid="A11",
        probe_path=test_files["probe_path"],
        frcmod=test_files["frcmod_path"],
        box_path=test_files["box_path"],
        size=10.0,
        ssbonds=[],
        at="gaff"
    )

def test_tleap_initialization():
    """TLeapクラスの初期化テスト"""
    tleap = TLeap()
    assert tleap.exe == os.getenv("TLEAP", "tleap"), "実行ファイルパスが正しくありません"
    assert not tleap.debug, "デバッグモードが正しく設定されていません"

    tleap_debug = TLeap(debug=True)
    assert tleap_debug.debug, "デバッグモードが有効化されていません"

def test_tleap_set_parameters(tleap: TLeap, test_files: dict):
    """パラメータ設定のテスト"""
    tleap_instance = tleap.set(
        cid="A11",
        probe_path=test_files["probe_path"],
        frcmod=test_files["frcmod_path"],
        box_path=test_files["box_path"],
        size=10.0,
        ssbonds=[],
        at="gaff"
    )

    assert tleap_instance.cid == "A11", "CIDが正しく設定されていません"
    assert tleap_instance.probe_path == test_files["probe_path"], "プローブパスが正しく設定されていません"
    assert tleap_instance.frcmod == test_files["frcmod_path"], "FRCMODパスが正しく設定されていません"
    assert tleap_instance.box_path == test_files["box_path"], "ボックスパスが正しく設定されていません"
    assert tleap_instance.size == 10.0, "サイズが正しく設定されていません"
    assert tleap_instance.ssbonds == [], "SSボンドが正しく設定されていません"
    assert tleap_instance.at == "gaff", "力場タイプが正しく設定されていません"

def test_tleap_run_normal(tleap_with_params: TLeap, tmp_path: Path):
    """正常系の実行テスト"""
    output_prefix = str(tmp_path / "test_output")
    
    try:
        result = tleap_with_params.run(output_prefix)
        assert result is not None, "実行結果がNoneです"
        assert hasattr(result, '_final_charge_value'), "final_charge_valueが存在しません"
        assert isinstance(result._final_charge_value, int), "final_charge_valueが整数ではありません"
    except Exception as e:
        pytest.fail(f"tleapの実行に失敗しました: {str(e)}")

def test_tleap_run_with_invalid_box(tleap_with_params: TLeap, tmp_path: Path):
    """異常系の実行テスト（無効なボックスファイル）"""
    tleap_with_params.box_path = Path("nonexistent.pdb")
    output_prefix = str(tmp_path / "test_output")
    
    with pytest.raises(subprocess.CalledProcessError) as exc_info:
        tleap_with_params.run(output_prefix)
    
    assert exc_info.value.returncode == 31, f"予期しないエラーコード: {exc_info.value.returncode}"

def test_tleap_run_with_ssbonds(tleap: TLeap, test_files: dict, tmp_path: Path):
    """SSボンド設定ありの実行テスト"""
    ssbonds: List[str] = []
    
    tleap_instance = tleap.set(
        cid="A11",
        probe_path=test_files["probe_path"],
        frcmod=test_files["frcmod_path"],
        box_path=test_files["box_path"],
        size=10.0,
        ssbonds=ssbonds,
        at="gaff"
    )
    
    assert tleap_instance.ssbonds == ssbonds, "SSボンドが正しく設定されていません"
    
    output_prefix = str(tmp_path / "test_output")
    try:
        result = tleap_instance.run(output_prefix)
        assert result is not None, "実行結果がNoneです"
    except Exception as e:
        pytest.fail(f"SSボンド設定なしのtleap実行に失敗しました: {str(e)}")
