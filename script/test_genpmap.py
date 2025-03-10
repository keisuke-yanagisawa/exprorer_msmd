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

# Basic grid-related fixtures
@pytest.fixture
def zero_grid_fixture():
    """Prepare basic grid data"""
    grid = gridData.Grid()
    grid.grid = np.zeros((2, 2, 2), dtype=np.float64)
    grid.origin = np.array([0.0, 0.0, 0.0])
    grid.delta = np.array([1.0, 1.0, 1.0])
    return grid

@pytest.fixture
def count_grid_fixture():
    """Prepare grid data for probability conversion test, probability values do not include 0"""
    grid = gridData.Grid()
    grid.grid = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]], dtype=np.float64)
    grid.origin = np.array([0.0, 0.0, 0.0])
    grid.delta = np.array([1.0, 1.0, 1.0])
    # Create a copy of grid data to prevent modifications
    grid.grid = grid.grid.copy()
    return grid

# Mask-related fixtures
@pytest.fixture
def mock_ref_struct_path():
    """Reference structure path for mock test
    Note: Actual file not needed, used only as a path object"""
    return Path("test_data/ref.pdb")

@pytest.fixture
def mask_fixture():
    """Prepare mask data"""
    return np.array([[[True, False], [True, False]], [[False, True], [False, True]]], dtype=bool)

# File I/O related fixtures
@pytest.fixture
def pmap_test_data(tmp_path):
    """Prepare data for PMAP conversion test"""
    grid_path = tmp_path / "test_grid.dx"
    ref_struct = Path("test_data/ref.pdb")
    
    # Create test grid file
    grid = gridData.Grid()
    grid.grid = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]], dtype=np.float64)
    grid.origin = np.array([0.0, 0.0, 0.0])
    grid.delta = np.array([1.0, 1.0, 1.0])
    grid.export(str(grid_path))
    
    return grid_path, ref_struct

# Test functions
@patch('script.utilities.GridUtil.gen_distance_grid')
def test_mask_generator_with_distance(mock_gen_distance, zero_grid_fixture, mock_ref_struct_path):
    """Test mask generation with distance threshold"""
    # Mock setup
    mock_grid = gridData.Grid()
    # Create grid with distance values (including boundary value test with 3.0)
    mock_grid.grid = np.array([
        [[2.0, 4.0], [3.0, 3.5]],
        [[2.9, 2.5], [3.0, 1.5]]
    ], dtype=np.float64)
    mock_grid.origin = np.array([0.0, 0.0, 0.0])
    mock_grid.delta = np.array([1.0, 1.0, 1.0])
    mock_gen_distance.return_value = mock_grid

    # Test execution (distance threshold 3.0)
    result = mask_generator(mock_ref_struct_path, zero_grid_fixture, distance=3.0)

    # Verify mock was called correctly
    mock_gen_distance.assert_called_once_with(zero_grid_fixture, mock_ref_struct_path)

    # Verify results
    # In mask_generator, mask.grid = mask.grid < distance, so
    # all values >= 3.0 become False
    expected_mask = np.array([
        [[True, False], [False, False]],  # 3.0 becomes False
        [[True, True], [False, True]]
    ], dtype=bool)
    np.testing.assert_array_equal(result.grid, expected_mask)

@patch('script.utilities.GridUtil.gen_distance_grid')
def test_mask_generator_without_distance(mock_gen_distance, zero_grid_fixture, mock_ref_struct_path):
    """Test mask generation without distance threshold"""
    # Mock setup
    mock_grid = gridData.Grid()
    mock_grid.grid = np.array([
        [[1.0, 2.0], [3.0, 4.0]],
        [[5.0, 6.0], [7.0, 8.0]]
    ])
    mock_grid.origin = np.array([0.0, 0.0, 0.0])
    mock_grid.delta = np.array([1.0, 1.0, 1.0])
    mock_gen_distance.return_value = mock_grid

    # Test execution
    result = mask_generator(mock_ref_struct_path, zero_grid_fixture)
    
    # Verify results
    # When distance=None, all finite values are masked
    expected_mask = np.array([
        [[True, True], [True, True]],
        [[True, True], [True, True]]
    ], dtype=bool)
    np.testing.assert_array_equal(result.grid, expected_mask)

@patch('script.utilities.GridUtil.gen_distance_grid')
def test_mask_generator_with_infinite_values(mock_gen_distance, zero_grid_fixture, mock_ref_struct_path):
    """Verify mask generation behavior for grid containing np.inf"""
    # Mock setup
    mock_grid = gridData.Grid()
    mock_grid.grid = np.array([
        [[1.0, np.inf], [2.0, np.inf]],
        [[np.inf, 3.0], [np.inf, 4.0]]
    ])
    mock_grid.origin = np.array([0.0, 0.0, 0.0])
    mock_grid.delta = np.array([1.0, 1.0, 1.0])
    mock_gen_distance.return_value = mock_grid

    # Test execution
    result = mask_generator(mock_ref_struct_path, zero_grid_fixture)
    
    # Verify results
    # Finite values become True, infinite values become False
    expected_mask = np.array([
        [[True, False], [True, False]],
        [[False, True], [False, True]]
    ], dtype=bool)
    np.testing.assert_array_equal(result.grid, expected_mask)


def test_convert_to_proba_with_mask_snapshot(count_grid_fixture, mask_fixture):
    """Test probability conversion with mask in snapshot mode"""
    result = convert_to_proba(count_grid_fixture, mask_fixture, normalize="snapshot", frames=2)
    # Verify masked values are divided by 2
    masked_values = result.grid[mask_fixture]
    np.testing.assert_array_almost_equal(masked_values, np.array([1, 3, 6, 8]) / 2)

def test_convert_to_proba_with_mask_total(count_grid_fixture, mask_fixture):
    """Test probability conversion with mask in total mode"""
    result = convert_to_proba(count_grid_fixture, mask_fixture, normalize="total")
    # Verify masked values sum to 1
    masked_values = result.grid[mask_fixture]
    assert np.abs(np.sum(masked_values) - 1.0) < 1e-10

def test_convert_to_proba_without_mask(count_grid_fixture):
    """Test probability conversion without mask"""
    result = convert_to_proba(count_grid_fixture)
    # Verify all values sum to 1
    assert 1.0 == pytest.approx(np.sum(result.grid), abs=1e-6)

def test_convert_to_proba_with_zero_values(count_grid_fixture, mask_fixture):
    """Test conversion of data containing zero values"""
    count_grid_fixture.grid[0, 0, 0] = 0.0
    result = convert_to_proba(count_grid_fixture, mask_fixture, normalize="total")
    # Verify masked values sum to 1
    masked_values = result.grid[mask_fixture]
    assert 1.0 == pytest.approx(np.sum(masked_values), abs=1e-6)
    # Verify unmasked regions are -1
    assert np.all(result.grid[~mask_fixture] == -1)
    # Verify zero values are replaced with 1e-10
    assert 1e-10 == pytest.approx(np.sum(result.grid[count_grid_fixture.grid == 0.0]), abs=1e-10)

def test_convert_to_gfe(tmp_path):
    """Test GFE conversion"""
    grid_path = tmp_path / "test_grid.dx"
    
    # Create test grid file
    pmap_grid = gridData.Grid()
    pmap_grid.grid = np.array([[[0.1, 0.2], [0.3, 0.4]], [[0.5, 0.6], [0.7, 0.8]]], dtype=np.float64)
    pmap_grid.origin = np.array([0.0, 0.0, 0.0])
    pmap_grid.delta = np.array([1.0, 1.0, 1.0])
    pmap_grid.export(str(grid_path))

    mean_proba = np.mean(pmap_grid.grid)
    temperature = 300
    
    # Execute GFE conversion
    gfe_path = convert_to_gfe(grid_path, mean_proba, temperature)
    
    # Verify results
    gfe_grid = gridData.Grid(gfe_path)
    RT = (constants.R / constants.calorie / constants.kilo) * temperature
    
    # Check value range (≤ 3 kcal/mol)
    assert np.all(gfe_grid.grid <= 3)
    
    # Verify negative GFE file is also generated
    inv_gfe_path = Path(gfe_path).parent / f"InvGFE_{Path(grid_path).name}"
    assert inv_gfe_path.exists()

@patch('script.genpmap.mask_generator')
@patch('script.genpmap.convert_to_proba')
def test_convert_to_pmap_snapshot(mock_convert_to_proba, mock_mask_generator, pmap_test_data):
    """Test PMAP conversion in snapshot mode"""
    grid_path, ref_struct = pmap_test_data
    
    # Mock setup
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

    # Test execution
    pmap_path = convert_to_pmap(grid_path, ref_struct, 3.0, normalize="snapshot", frames=10)

    # Verify mocks were called correctly
    mock_mask_generator.assert_called_once()
    mock_convert_to_proba.assert_called_once()

    # Verify results
    assert os.path.exists(pmap_path)
    result_grid = gridData.Grid(pmap_path)
    assert result_grid.grid.shape == (2, 2, 2)

@patch('script.genpmap.mask_generator')
@patch('script.genpmap.convert_to_proba')
def test_convert_to_pmap_total(mock_convert_to_proba, mock_mask_generator, pmap_test_data):
    """Test PMAP conversion in total mode"""
    grid_path, ref_struct = pmap_test_data
    
    # Mock setup
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

    # Test execution
    pmap_path = convert_to_pmap(grid_path, ref_struct, 3.0, normalize="total")

    # Verify results
    result_grid = gridData.Grid(pmap_path)
    assert np.abs(np.sum(result_grid.grid[mock_mask.grid]) - 1.0) < 1e-10

def test_parse_basic_range():
    """Test parsing basic range specification"""
    start, stop, offset = parse_snapshot_setting("1-100")
    assert start == "1"
    assert stop == "100"
    assert offset == "1"  # Default value

def test_parse_with_offset():
    """Test parsing range specification with offset"""
    start, stop, offset = parse_snapshot_setting("1-100:2")
    assert start == "1"
    assert stop == "100"
    assert offset == "2"

def test_parse_single_frame():
    """Test parsing single frame specification"""
    start, stop, offset = parse_snapshot_setting("1-1")
    assert start == "1"
    assert stop == "1"
    assert offset == "1"

def test_parse_invalid_format():
    """Test invalid format input"""
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
    """Prepare test data for PMAP generation"""
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
    """Test basic PMAP generation"""
    setting_general, setting_input, setting_pmap, traj, top = gen_pmap_test_data
    
    # Mock setup
    mock_get_attr.return_value = np.array([[1.0, 1.0, 1.0]])
    mock_cpptraj_instance = MagicMock()
    mock_cpptraj_instance.maps = [{"grid": Path("test_data/map1.dx")}, {"grid": Path("test_data/map2.dx")}]
    mock_cpptraj_instance.frames = 50
    mock_cpptraj_instance.last_volume = 1000.0
    mock_cpptraj.return_value = mock_cpptraj_instance
    mock_convert_to_pmap.side_effect = ["test_data/pmap1.dx", "test_data/pmap2.dx"]

    # Test execution
    pmap_paths = gen_pmap(
        tmp_path,
        setting_general,
        setting_input,
        setting_pmap,
        traj,
        top
    )

    # Verify results
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
    """Test PMAP generation in GFE mode"""
    setting_general, setting_input, setting_pmap, traj, top = gen_pmap_test_data
    
    # Change settings to GFE mode
    setting_pmap["normalization"] = "GFE"
    
    # Mock setup
    mock_get_attr.return_value = np.array([[1.0, 1.0, 1.0]])
    mock_cpptraj_instance = MagicMock()
    mock_cpptraj_instance.maps = [{"grid": Path("test_data/map1.dx"), "num_probe_atoms": 100}]
    mock_cpptraj_instance.frames = 50
    mock_cpptraj_instance.last_volume = 1000.0
    mock_cpptraj.return_value = mock_cpptraj_instance
    mock_estimate_volume.return_value = 100.0
    mock_convert_to_pmap.return_value = "test_data/pmap1.dx"
    mock_convert_to_gfe.return_value = "test_data/gfe1.dx"

    # Test execution
    pmap_paths = gen_pmap(
        tmp_path,
        setting_general,
        setting_input,
        setting_pmap,
        traj,
        top
    )

    # Verify results
    assert len(pmap_paths) == 1
    mock_convert_to_gfe.assert_called_once()
    # Verify mean_proba calculation is correct
    expected_mean_proba = 100 / (1000.0 - 100.0)
    mock_convert_to_gfe.assert_called_with("test_data/pmap1.dx", expected_mean_proba, temperature=300)
