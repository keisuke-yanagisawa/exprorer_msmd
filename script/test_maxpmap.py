import os
from pathlib import Path
import pytest
import numpy as np
from gridData import Grid

from script.maxpmap import grid_max, gen_max_pmap

# テストデータのパス定義
TEST_DATA_DIR = Path("script/test_data")
PMAP1_PATH = TEST_DATA_DIR / "pmap1.dx"
PMAP2_PATH = TEST_DATA_DIR / "pmap2.dx"
SMALL_GRID_PATH = TEST_DATA_DIR / "small_grid.dx"

@pytest.fixture
def grid1():
    """pmap1.dxからグリッドを読み込む"""
    return Grid(str(PMAP1_PATH))

@pytest.fixture
def grid2():
    """pmap2.dxからグリッドを読み込む"""
    return Grid(str(PMAP2_PATH))

@pytest.fixture
def small_grid():
    """small_grid.dxからグリッドを読み込む"""
    return Grid(str(SMALL_GRID_PATH))


class TestGridMax:
    """grid_max関数のテストクラス"""
    
    def test_maxpmap(self, grid1, grid2):
        """基本的な最大値の計算テスト"""
        maxpmap = grid_max([grid1, grid2])
        assert grid1.grid[1, 5, 2] == pytest.approx(maxpmap.grid[1, 5, 2])  # -1
        assert grid1.grid[38, 59, 27] == pytest.approx(maxpmap.grid[38, 59, 27])  # 0
        assert grid1.grid[44, 20, 32] == pytest.approx(
            maxpmap.grid[44, 20, 32]
        )  # not -1 & not 0 for grid1, grid2[44,20,32] = 0

    def test_empty_grid_list(self):
        """空のグリッドリストに対するテスト"""
        with pytest.raises(ValueError, match="Empty grid list"):
            grid_max([])

    def test_single_grid(self, small_grid):
        """1つのグリッドに対するテスト"""
        result = grid_max([small_grid])
        assert result.grid == pytest.approx(small_grid.grid)

    def test_different_size_grids(self, grid1, small_grid):
        """異なるサイズのグリッドに対するテスト"""
        with pytest.raises(ValueError, match="Grids have different sizes"):
            grid_max([grid1, small_grid])

    def test_different_origin_grids(self, small_grid):
        """異なる原点を持つグリッドに対するテスト"""
        grid_modified = Grid(str(SMALL_GRID_PATH))
        grid_modified.origin = np.array([1.0, 1.0, 1.0])
        
        with pytest.raises(ValueError, match="Grids have different origins"):
            grid_max([small_grid, grid_modified])

    def test_different_delta_grids(self, small_grid):
        """異なるデルタ値を持つグリッドに対するテスト"""
        grid_modified = Grid(str(SMALL_GRID_PATH))
        grid_modified.delta = np.array([2.0, 2.0, 2.0])
        
        with pytest.raises(ValueError, match="Grids have different deltas"):
            grid_max([small_grid, grid_modified])


class TestGenMaxPmap:
    """gen_max_pmap関数のテストクラス"""
    
    def test_gen_max_pmap(self, small_grid, tmp_path):
        """gen_max_pmap関数の基本機能テスト"""
        temp_output = tmp_path / "temp_max.dx"
        input_paths = [str(SMALL_GRID_PATH)]
        output_path = gen_max_pmap(input_paths, str(temp_output))
        
        assert Path(output_path).exists()
        result_grid = Grid(output_path)
        original_grid = Grid(input_paths[0])
        assert result_grid.grid == pytest.approx(original_grid.grid)

    def test_empty_input(self):
        """空の入力リストに対するgen_max_pmapのテスト"""
        with pytest.raises(ValueError, match="No input files provided"):
            gen_max_pmap([], "output.dx")