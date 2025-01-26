import pytest
import numpy as np
import gridData
from pathlib import Path
import os
from scipy import constants
from unittest.mock import patch, MagicMock
import numpy.testing as npt

from script.genpmap import (
    mask_generator,
    convert_to_proba,
    convert_to_gfe,
    convert_to_pmap,
    parse_snapshot_setting,
    gen_pmap
)

# 基本的なグリッド関連のfixture
@pytest.fixture
def zero_grid_fixture():
    """基本的なグリッドデータの準備"""
    grid = gridData.Grid()
    grid.grid = np.zeros((2, 2, 2), dtype=np.float64)
    grid.origin = np.array([0.0, 0.0, 0.0])
    grid.delta = np.array([1.0, 1.0, 1.0])
    return grid

@pytest.fixture
def count_grid_fixture():
    """確率変換テスト用のグリッドデータ準備、確率値は0を含まない"""
    grid = gridData.Grid()
    grid.grid = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]], dtype=np.float64)
    grid.origin = np.array([0.0, 0.0, 0.0])
    grid.delta = np.array([1.0, 1.0, 1.0])
    # グリッドデータのコピーを作成して変更を防ぐ
    grid.grid = grid.grid.copy()
    return grid

# マスク関連のfixture
@pytest.fixture
def mock_ref_struct_path():
    """モックテスト用の参照構造パス
    Note: 実際のファイルは不要で、パスオブジェクトとしてのみ使用"""
    return Path("test_data/ref.pdb")

@pytest.fixture
def mask_fixture():
    """マスクデータの準備"""
    return np.array([[[True, False], [True, False]], [[False, True], [False, True]]], dtype=bool)

# ファイルI/O関連のfixture
@pytest.fixture
def pmap_test_data(tmp_path):
    """PMAP変換テスト用のデータ準備"""
    grid_path = tmp_path / "test_grid.dx"
    ref_struct = Path("test_data/ref.pdb")
    
    # テスト用のグリッドファイルを作成
    grid = gridData.Grid()
    grid.grid = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]], dtype=np.float64)
    grid.origin = np.array([0.0, 0.0, 0.0])
    grid.delta = np.array([1.0, 1.0, 1.0])
    grid.export(str(grid_path))
    
    return grid_path, ref_struct

# テスト関数
@patch('script.utilities.GridUtil.gen_distance_grid')
def test_mask_generator_with_distance(mock_gen_distance, zero_grid_fixture, mock_ref_struct_path):
    """距離閾値ありのマスク生成テスト"""
    # モックの設定
    mock_grid = gridData.Grid()
    # 距離値を持つグリッドを作成（3.0を含む境界値テスト）
    mock_grid.grid = np.array([
        [[2.0, 4.0], [3.0, 3.5]],
        [[2.9, 2.5], [3.0, 1.5]]
    ], dtype=np.float64)
    mock_grid.origin = np.array([0.0, 0.0, 0.0])
    mock_grid.delta = np.array([1.0, 1.0, 1.0])
    mock_gen_distance.return_value = mock_grid

    # テスト実行（距離閾値3.0）
    result = mask_generator(mock_ref_struct_path, zero_grid_fixture, distance=3.0)

    # モックが正しく呼び出されたことを確認
    mock_gen_distance.assert_called_once_with(zero_grid_fixture, mock_ref_struct_path)

    # 結果の検証
    # mask_generator内で mask.grid = mask.grid < distance としているため、
    # 3.0以上の値はすべてFalseになる
    expected_mask = np.array([
        [[True, False], [False, False]],  # 3.0はFalseになる
        [[True, True], [False, True]]
    ], dtype=bool)
    np.testing.assert_array_equal(result.grid, expected_mask)
@patch('script.utilities.GridUtil.gen_distance_grid')
def test_mask_generator_without_distance(mock_gen_distance, zero_grid_fixture, mock_ref_struct_path):
    """距離閾値なしのマスク生成テスト"""
    # モックの設定
    mock_grid = gridData.Grid()
    mock_grid.grid = np.array([
        [[1.0, 2.0], [3.0, 4.0]],
        [[5.0, 6.0], [7.0, 8.0]]
    ])
    mock_grid.origin = np.array([0.0, 0.0, 0.0])
    mock_grid.delta = np.array([1.0, 1.0, 1.0])
    mock_gen_distance.return_value = mock_grid

    # テスト実行
    result = mask_generator(mock_ref_struct_path, zero_grid_fixture)
    
    # 結果の検証
    # distance=Noneの場合、全ての有限値がマスクされる
    expected_mask = np.array([
        [[True, True], [True, True]],
        [[True, True], [True, True]]
    ], dtype=bool)
    np.testing.assert_array_equal(result.grid, expected_mask)

@patch('script.utilities.GridUtil.gen_distance_grid')
def test_mask_generator_with_infinite_values(mock_gen_distance, zero_grid_fixture, mock_ref_struct_path):
    """np.infを含むグリッドに対するマスク生成の挙動を確認"""
    # モックの設定
    mock_grid = gridData.Grid()
    mock_grid.grid = np.array([
        [[1.0, np.inf], [2.0, np.inf]],
        [[np.inf, 3.0], [np.inf, 4.0]]
    ])
    mock_grid.origin = np.array([0.0, 0.0, 0.0])
    mock_grid.delta = np.array([1.0, 1.0, 1.0])
    mock_gen_distance.return_value = mock_grid

    # テスト実行
    result = mask_generator(mock_ref_struct_path, zero_grid_fixture)
    
    # 結果の検証
    # 有限の値はTrue、無限大の値はFalseになる
    expected_mask = np.array([
        [[True, False], [True, False]],
        [[False, True], [False, True]]
    ], dtype=bool)
    np.testing.assert_array_equal(result.grid, expected_mask)


def test_convert_to_proba_with_mask_snapshot(count_grid_fixture, mask_fixture):
    """マスクありでsnapshotモードのテスト"""
    result = convert_to_proba(count_grid_fixture, mask_fixture, normalize="snapshot", frames=2)
    # マスクされた領域の値が2で割られることを確認
    masked_values = result.grid[mask_fixture]
    np.testing.assert_array_almost_equal(masked_values, np.array([1, 3, 6, 8]) / 2)

def test_convert_to_proba_with_mask_total(count_grid_fixture, mask_fixture):
    """マスクありでtotalモードのテスト"""
    result = convert_to_proba(count_grid_fixture, mask_fixture, normalize="total")
    # マスクされた領域の値が合計で1になることを確認
    masked_values = result.grid[mask_fixture]
    assert np.abs(np.sum(masked_values) - 1.0) < 1e-10

def test_convert_to_proba_without_mask(count_grid_fixture):
    """マスクなしのテスト"""
    result = convert_to_proba(count_grid_fixture)
    # 全体の値が合計で1になることを確認
    assert 1.0 == pytest.approx(np.sum(result.grid), abs=1e-6)

def test_convert_to_proba_with_mask_gfe(count_grid_fixture, mask_fixture):
    """マスクありでGFEモードのテスト、result自体はsnapshotと同じ結果が出力される。"""
    result = convert_to_proba(count_grid_fixture, mask_fixture, normalize="GFE", frames=2)
    # マスクされた領域の値がフレーム数 (=2) で割られることを確認（snapshotと同じ動作）
    masked_values = result.grid[mask_fixture]
    assert np.array([1, 3, 6, 8]) / 2 == pytest.approx(masked_values, abs=1e-6)

def test_convert_to_proba_with_zero_values(count_grid_fixture, mask_fixture):
    """ゼロ値を含むデータの変換テスト"""
    count_grid_fixture.grid[0, 0, 0] = 0.0
    result = convert_to_proba(count_grid_fixture, mask_fixture, normalize="total")
    # マスクされた領域の値が合計で1になることを確認
    masked_values = result.grid[mask_fixture]
    assert 1.0 == pytest.approx(np.sum(masked_values), abs=1e-6)
    # マスクされていない領域が-1になることを確認
    assert np.all(result.grid[~mask_fixture] == -1)
    # ゼロ値が1e-10に置換されていることを確認
    assert 1e-10 == pytest.approx(np.sum(result.grid[count_grid_fixture.grid == 0.0]), abs=1e-10)

def test_convert_to_gfe(tmp_path):
    """GFE変換のテスト"""
    grid_path = tmp_path / "test_grid.dx"
    
    # テスト用のグリッドファイルを作成
    grid = gridData.Grid()
    grid.grid = np.array([[[0.1, 0.2], [0.3, 0.4]], [[0.5, 0.6], [0.7, 0.8]]], dtype=np.float64)
    grid.origin = np.array([0.0, 0.0, 0.0])
    grid.delta = np.array([1.0, 1.0, 1.0])
    grid.export(str(grid_path))

    mean_proba = 0.5
    temperature = 300
    
    # GFE変換を実行
    gfe_path = convert_to_gfe(grid_path, mean_proba, temperature)
    
    # 結果の検証
    gfe_grid = gridData.Grid(gfe_path)
    RT = (constants.R / constants.calorie / constants.kilo) * temperature
    
    # 値の範囲チェック（3 kcal/mol以下）
    assert np.all(gfe_grid.grid <= 3)
    
    # 負のGFEファイルも生成されることを確認
    inv_gfe_path = Path(gfe_path).parent / f"InvGFE_{Path(grid_path).name}"
    assert inv_gfe_path.exists()

@patch('script.genpmap.mask_generator')
@patch('script.genpmap.convert_to_proba')
def test_convert_to_pmap_snapshot(mock_convert_to_proba, mock_mask_generator, pmap_test_data):
    """snapshotモードでのPMAP変換テスト"""
    grid_path, ref_struct = pmap_test_data
    
    # モックの設定
    mock_mask = gridData.Grid()
    mock_mask.grid = np.array([[[True, False], [True, False]], [[False, True], [False, True]]], dtype=bool)
    mock_mask.origin = np.array([0.0, 0.0, 0.0])
    mock_mask.delta = np.array([1.0, 1.0, 1.0])
    mock_mask_generator.return_value = mock_mask
    
    mock_pmap = gridData.Grid()
    mock_pmap.grid = np.array([[[0.1, 0.0], [0.2, 0.0]], [[0.0, 0.3], [0.0, 0.4]]], dtype=np.float64)
    mock_pmap.origin = np.array([0.0, 0.0, 0.0])
    mock_pmap.delta = np.array([1.0, 1.0, 1.0])
    mock_convert_to_proba.return_value = mock_pmap

    # テスト実行
    pmap_path = convert_to_pmap(grid_path, ref_struct, 3.0, normalize="snapshot", frames=10)

    # モックが正しく呼び出されたことを確認
    mock_mask_generator.assert_called_once()
    mock_convert_to_proba.assert_called_once()

    # 結果の検証
    assert os.path.exists(pmap_path)
    result_grid = gridData.Grid(pmap_path)
    assert result_grid.grid.shape == (2, 2, 2)

@patch('script.genpmap.mask_generator')
@patch('script.genpmap.convert_to_proba')
def test_convert_to_pmap_total(mock_convert_to_proba, mock_mask_generator, pmap_test_data):
    """totalモードでのPMAP変換テスト"""
    grid_path, ref_struct = pmap_test_data
    
    # モックの設定
    mock_mask = gridData.Grid()
    mock_mask.grid = np.array([[[True, False], [True, False]], [[False, True], [False, True]]], dtype=bool)
    mock_mask.origin = np.array([0.0, 0.0, 0.0])
    mock_mask.delta = np.array([1.0, 1.0, 1.0])
    mock_mask_generator.return_value = mock_mask
    
    mock_pmap = gridData.Grid()
    mock_pmap.grid = np.array([[[0.25, 0.0], [0.25, 0.0]], [[0.0, 0.25], [0.0, 0.25]]], dtype=np.float64)
    mock_pmap.origin = np.array([0.0, 0.0, 0.0])
    mock_pmap.delta = np.array([1.0, 1.0, 1.0])
    mock_convert_to_proba.return_value = mock_pmap

    # テスト実行
    pmap_path = convert_to_pmap(grid_path, ref_struct, 3.0, normalize="total")

    # 結果の検証
    result_grid = gridData.Grid(pmap_path)
    assert np.abs(np.sum(result_grid.grid[mock_mask.grid]) - 1.0) < 1e-10

def test_parse_basic_range():
    """基本的な範囲指定のパース"""
    start, stop, offset = parse_snapshot_setting("1-100")
    assert start == "1"
    assert stop == "100"
    assert offset == "1"  # デフォルト値

def test_parse_with_offset():
    """オフセット付きの範囲指定のパース"""
    start, stop, offset = parse_snapshot_setting("1-100:2")
    assert start == "1"
    assert stop == "100"
    assert offset == "2"

def test_parse_single_frame():
    """単一フレームの指定のパース"""
    start, stop, offset = parse_snapshot_setting("1-1")
    assert start == "1"
    assert stop == "1"
    assert offset == "1"

def test_parse_invalid_format():
    """無効な形式の入力テスト"""
    with pytest.raises(ValueError):
        parse_snapshot_setting("invalid")
    with pytest.raises(ValueError):
        parse_snapshot_setting("1:2:3")
    with pytest.raises(ValueError):
        parse_snapshot_setting("1-")
    with pytest.raises(ValueError):
        parse_snapshot_setting("-1")
    with pytest.raises(ValueError):
        parse_snapshot_setting(":")

@pytest.fixture
def gen_pmap_test_data(tmp_path):
    """PMAP生成テスト用のデータ準備"""
    setting_general = {
        "name": "test_simulation"
    }
    setting_input = {
        "protein": {
            "pdb": "test_data/ref.pdb"
        },
        "probe": {
            "cid": "TEST"
        }
    }
    setting_pmap = {
        "snapshot": "1-100:2",
        "maps": ["density", "occupancy"],
        "map_size": 20,
        "valid_dist": 3.0,
        "normalization": "snapshot"
    }
    traj = Path("test_data/traj.xtc")
    top = Path("test_data/top.pdb")
    
    return setting_general, setting_input, setting_pmap, traj, top

@patch('script.genpmap.uPDB.get_structure')
@patch('script.genpmap.uPDB.get_attr')
@patch('script.genpmap.Cpptraj')
@patch('script.genpmap.convert_to_pmap')
def test_gen_pmap_basic(mock_convert_to_pmap, mock_cpptraj, mock_get_attr, mock_get_structure, tmp_path, gen_pmap_test_data):
    """基本的なPMAP生成テスト"""
    setting_general, setting_input, setting_pmap, traj, top = gen_pmap_test_data
    
    # モックの設定
    mock_get_attr.return_value = np.array([[1.0, 1.0, 1.0]])
    mock_cpptraj_instance = MagicMock()
    mock_cpptraj_instance.maps = [{"grid": Path("test_data/map1.dx")}, {"grid": Path("test_data/map2.dx")}]
    mock_cpptraj_instance.frames = 50
    mock_cpptraj_instance.last_volume = 1000.0
    mock_cpptraj.return_value = mock_cpptraj_instance
    mock_convert_to_pmap.side_effect = ["test_data/pmap1.dx", "test_data/pmap2.dx"]

    # テスト実行
    pmap_paths = gen_pmap(
        tmp_path,
        setting_general,
        setting_input,
        setting_pmap,
        traj,
        top
    )

    # 結果の検証
    assert len(pmap_paths) == 2
    mock_cpptraj_instance.set.assert_called_once()
    mock_cpptraj_instance.run.assert_called_once()
    assert mock_convert_to_pmap.call_count == 2

@patch('script.genpmap.uPDB.get_structure')
@patch('script.genpmap.uPDB.get_attr')
@patch('script.genpmap.uPDB.estimate_exclute_volume')
@patch('script.genpmap.Cpptraj')
@patch('script.genpmap.convert_to_pmap')
@patch('script.genpmap.convert_to_gfe')
def test_gen_pmap_gfe(mock_convert_to_gfe, mock_convert_to_pmap, mock_cpptraj, 
                     mock_estimate_volume, mock_get_attr, mock_get_structure, tmp_path, gen_pmap_test_data):
    """GFEモードでのPMAP生成テスト"""
    setting_general, setting_input, setting_pmap, traj, top = gen_pmap_test_data
    
    # 設定をGFEモードに変更
    setting_pmap["normalization"] = "GFE"
    
    # モックの設定
    mock_get_attr.return_value = np.array([[1.0, 1.0, 1.0]])
    mock_cpptraj_instance = MagicMock()
    mock_cpptraj_instance.maps = [{"grid": Path("test_data/map1.dx"), "num_probe_atoms": 100}]
    mock_cpptraj_instance.frames = 50
    mock_cpptraj_instance.last_volume = 1000.0
    mock_cpptraj.return_value = mock_cpptraj_instance
    mock_estimate_volume.return_value = 100.0
    mock_convert_to_pmap.return_value = "test_data/pmap1.dx"
    mock_convert_to_gfe.return_value = "test_data/gfe1.dx"

    # テスト実行
    pmap_paths = gen_pmap(
        tmp_path,
        setting_general,
        setting_input,
        setting_pmap,
        traj,
        top
    )

    # 結果の検証
    assert len(pmap_paths) == 1
    mock_convert_to_gfe.assert_called_once()
    # mean_probaの計算が正しいことを確認
    expected_mean_proba = 100 / (1000.0 - 100.0)
    mock_convert_to_gfe.assert_called_with("test_data/pmap1.dx", expected_mean_proba, temperature=300)
