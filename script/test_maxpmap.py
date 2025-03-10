import os
from pathlib import Path
import pytest
import numpy as np
from gridData import Grid

from script.maxpmap import grid_max, gen_max_pmap

# Define test data paths
TEST_DATA_DIR = Path("script/test_data")
PMAP1_PATH = TEST_DATA_DIR / "pmap1.dx"
PMAP2_PATH = TEST_DATA_DIR / "pmap2.dx"
SMALL_GRID_PATH = TEST_DATA_DIR / "small_grid.dx"

@pytest.fixture
def grid1():
    """Load grid from pmap1.dx"""
    return Grid(str(PMAP1_PATH))

@pytest.fixture
def grid2():
    """Load grid from pmap2.dx"""
    return Grid(str(PMAP2_PATH))

@pytest.fixture
def small_grid():
    """Load grid from small_grid.dx"""
    return Grid(str(SMALL_GRID_PATH))


class TestGridMax:
    """Test class for grid_max function"""
    
    def test_maxpmap(self, grid1, grid2):
        """Test basic maximum value calculation"""
        maxpmap = grid_max([grid1, grid2])
        assert grid1.grid[1, 5, 2] == pytest.approx(maxpmap.grid[1, 5, 2])  # -1
        assert grid1.grid[38, 59, 27] == pytest.approx(maxpmap.grid[38, 59, 27])  # 0
        assert grid1.grid[44, 20, 32] == pytest.approx(
            maxpmap.grid[44, 20, 32]
        )  # not -1 & not 0 for grid1, grid2[44,20,32] = 0

    def test_empty_grid_list(self):
        """Test for empty grid list"""
        with pytest.raises(ValueError, match="Empty grid list"):
            grid_max([])

    def test_single_grid(self, small_grid):
        """Test for single grid"""
        result = grid_max([small_grid])
        assert result.grid == pytest.approx(small_grid.grid)

    def test_different_size_grids(self, grid1, small_grid):
        """Test for grids with different sizes"""
        with pytest.raises(ValueError, match="Grids have different sizes"):
            grid_max([grid1, small_grid])

    def test_different_origin_grids(self, small_grid):
        """Test for grids with different origins"""
        grid_modified = Grid(str(SMALL_GRID_PATH))
        grid_modified.origin = np.array([1.0, 1.0, 1.0])
        
        with pytest.raises(ValueError, match="Grids have different origins"):
            grid_max([small_grid, grid_modified])

    def test_different_delta_grids(self, small_grid):
        """Test for grids with different delta values"""
        grid_modified = Grid(str(SMALL_GRID_PATH))
        grid_modified.delta = np.array([2.0, 2.0, 2.0])
        
        with pytest.raises(ValueError, match="Grids have different deltas"):
            grid_max([small_grid, grid_modified])


class TestGenMaxPmap:
    """Test class for gen_max_pmap function"""
    
    def test_gen_max_pmap(self, small_grid, tmp_path):
        """Test basic functionality of gen_max_pmap function"""
        temp_output = tmp_path / "temp_max.dx"
        input_paths = [str(SMALL_GRID_PATH)]
        output_path = gen_max_pmap(input_paths, str(temp_output))
        
        assert Path(output_path).exists()
        result_grid = Grid(output_path)
        original_grid = Grid(input_paths[0])
        assert result_grid.grid == pytest.approx(original_grid.grid)

    def test_empty_input(self):
        """Test gen_max_pmap with empty input list"""
        with pytest.raises(ValueError, match="No input files provided"):
            gen_max_pmap([], "output.dx")