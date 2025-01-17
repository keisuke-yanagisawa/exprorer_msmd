import os
from pathlib import Path
from unittest import TestCase

import gridData
import numpy as np

from script.genpmap import mask_generator
from script.utilities import GridUtil


class TestMaskGenerator(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestMaskGenerator, self).__init__(*args, **kwargs)
        __testdata_dir = Path("script/test_data")
        self.grid = gridData.Grid()
        self.grid.load(f"{__testdata_dir}/small_grid.dx")
        self.ref_struct = __testdata_dir / "tripeptide.pdb"

    def test_generate_mask_distance_none(self):
        mask = mask_generator(self.ref_struct, self.grid, distance=None)
        self.assertTrue(np.all(mask.grid))
        mask = GridUtil.gen_distance_grid(self.grid, self.ref_struct)
        # def mask_generator(ref_struct, reference_grid, distance=None):
        #     mask = GridUtil.gen_distance_grid(reference_grid, ref_struct)
        #     # print(np.max(mask.grid), np.min(mask.grid), distance)
        #     if distance is not None:
        #         mask.grid = mask.grid < distance
        #     else:
        #         mask.grid = mask.grid < np.inf
        #     return mask

    def test_generate_mask_distance_with_distance_threshold(self):
        mask = mask_generator(self.ref_struct, self.grid, distance=18)
        self.assertTrue(mask.grid[0, 0, 0])
        self.assertTrue(mask.grid[0, 0, 1])
        self.assertFalse(mask.grid[0, 1, 0])
        self.assertTrue(mask.grid[0, 1, 1])
        self.assertFalse(mask.grid[1, 0, 0])
        self.assertFalse(mask.grid[1, 0, 1])
        self.assertFalse(mask.grid[1, 1, 0])
        self.assertFalse(mask.grid[1, 1, 1])

    def test_convert_to_proba_snapshot_normalization(self):
        """Test probability conversion with snapshot normalization"""
        from script.genpmap import convert_to_proba
        test_grid = gridData.Grid()
        test_grid.grid = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
        mask = np.ones_like(test_grid.grid, dtype=bool)
        frames = 2
        
        result = convert_to_proba(test_grid, mask, normalize="snapshot", frames=frames)
        # Values should be divided by number of frames
        expected = np.array([[[0.5, 1.0], [1.5, 2.0]], [[2.5, 3.0], [3.5, 4.0]]])
        np.testing.assert_array_almost_equal(result.grid, expected)

    def test_convert_to_proba_GFE_normalization(self):
        """Test probability conversion with GFE normalization"""
        from script.genpmap import convert_to_proba
        test_grid = gridData.Grid()
        test_grid.grid = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
        mask = np.ones_like(test_grid.grid, dtype=bool)
        frames = 2
        
        result = convert_to_proba(test_grid, mask, normalize="GFE", frames=frames)
        # GFE normalization should behave like snapshot normalization
        expected = np.array([[[0.5, 1.0], [1.5, 2.0]], [[2.5, 3.0], [3.5, 4.0]]])
        np.testing.assert_array_almost_equal(result.grid, expected)

    def test_convert_to_proba_others_normalization(self):
        """Test probability conversion with total normalization"""
        from script.genpmap import convert_to_proba
        test_grid = gridData.Grid()
        test_grid.grid = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
        mask = np.ones_like(test_grid.grid, dtype=bool)
        
        result = convert_to_proba(test_grid, mask, normalize="total")
        # Values should be divided by sum of all values (36.0)
        expected = test_grid.grid / 36.0
        np.testing.assert_array_almost_equal(result.grid, expected)
        
        # Test with partial mask
        mask[0, 0, 0] = False
        result = convert_to_proba(test_grid, mask, normalize="total")
        # Sum should exclude masked value (sum = 35.0)
        masked_sum = np.sum(test_grid.grid[mask])
        expected = np.full_like(test_grid.grid, -1.0)  # -1 for masked regions
        expected[mask] = test_grid.grid[mask] / masked_sum
        np.testing.assert_array_almost_equal(result.grid, expected)

    def test_convert_to_gfe(self):
        """Test conversion of probability map to GFE map"""
        import tempfile
        from script.genpmap import convert_to_gfe
        
        # Create a temporary grid file
        test_grid = gridData.Grid()
        test_grid.grid = np.array([[[0.1, 0.2], [0.0, 0.4]], [[0.0001, 0.6], [0.7, 0.8]]])
        temp_file = tempfile.mkstemp(suffix='.dx')[1]
        test_grid.export(temp_file)
        
        # Test GFE conversion
        mean_proba = 0.3
        temperature = 300
        gfe_path = convert_to_gfe(temp_file, mean_proba, temperature)
        
        # Load and check GFE grid
        gfe_grid = gridData.Grid(gfe_path)
        # Check values are capped at 3.0
        self.assertTrue(np.all(gfe_grid.grid <= 3.0))
        # Check small probabilities are handled correctly
        self.assertAlmostEqual(gfe_grid.grid[1, 0, 0], 3.0)  # Very small probability (0.0001)
        
        # Load and check inverse GFE grid
        inv_gfe_path = os.path.dirname(gfe_path) + "/InvGFE_" + os.path.basename(temp_file)
        inv_gfe_grid = gridData.Grid(inv_gfe_path)
        np.testing.assert_array_almost_equal(inv_gfe_grid.grid, -gfe_grid.grid)
        
        # Cleanup
        os.remove(temp_file)
        os.remove(gfe_path)
        os.remove(inv_gfe_path)

    def test_parse_snapshot_setting(self):
        """Test parsing of snapshot settings string"""
        from script.genpmap import parse_snapshot_setting
        
        # Test basic range
        start, stop, offset = parse_snapshot_setting("1-10")
        self.assertEqual((start, stop, offset), ("1", "10", "1"))
        
        # Test range with offset
        start, stop, offset = parse_snapshot_setting("1-10:2")
        self.assertEqual((start, stop, offset), ("1", "10", "2"))
        
        # Test larger numbers
        start, stop, offset = parse_snapshot_setting("100-200:5")
        self.assertEqual((start, stop, offset), ("100", "200", "5"))

    def test_convert_to_pmap_snapshot_normalization(self):
        """Test conversion to probability map with snapshot normalization"""
        import tempfile
        from script.genpmap import convert_to_pmap
        
        # Create test grid
        test_grid = gridData.Grid()
        test_grid.grid = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
        temp_grid = tempfile.mkstemp(suffix='.dx')[1]
        test_grid.export(temp_grid)
        
        # Convert to pmap
        pmap_path = convert_to_pmap(
            Path(temp_grid),
            self.ref_struct,
            valid_distance=18.0,
            normalize="snapshot",
            frames=2
        )
        
        # Load and check pmap
        pmap = gridData.Grid(pmap_path)
        # Check some values are properly normalized
        self.assertAlmostEqual(np.max(pmap.grid[pmap.grid > -1]), 4.0)  # Max value / frames
        
        # Cleanup
        os.remove(temp_grid)
        os.remove(pmap_path)

    def test_convert_to_pmap_others_normalization(self):
        """Test conversion to probability map with total normalization"""
        import tempfile
        from script.genpmap import convert_to_pmap
        
        # Create test grid
        test_grid = gridData.Grid()
        test_grid.grid = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]])
        temp_grid = tempfile.mkstemp(suffix='.dx')[1]
        test_grid.export(temp_grid)
        
        # Convert to pmap with total normalization
        pmap_path = convert_to_pmap(
            Path(temp_grid),
            self.ref_struct,
            valid_distance=18.0,
            normalize="total"
        )
        
        # Load and check pmap
        pmap = gridData.Grid(pmap_path)
        valid_values = pmap.grid[pmap.grid > -1]
        self.assertAlmostEqual(np.sum(valid_values), 1.0)  # Sum should be 1
        
        # Cleanup
        os.remove(temp_grid)
        os.remove(pmap_path)

    def test_gen_pmap_snapshot_normalization(self):
        """Test full pmap generation with snapshot normalization"""
        from script.genpmap import gen_pmap
        import tempfile
        
        # Create minimal settings
        setting_general = {"name": "test"}
        setting_input = {
            "protein": {"pdb": str(self.ref_struct)},
            "probe": {"cid": "TEST"}
        }
        setting_pmap = {
            "snapshot": "1-10",
            "valid_dist": 18.0,
            "map_size": 80,
            "normalization": "snapshot",
            "maps": [{"suffix": "test", "selector": "!@H*"}]
        }
        
        # Test pmap generation
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = gen_pmap(
                Path(tmpdir),
                setting_general,
                setting_input,
                setting_pmap,
                Path(tmpdir) / "traj.nc",
                Path(tmpdir) / "top.parm7"
            )
            self.assertTrue(len(paths) > 0)

    def test_gen_pmap_gfe_normalization(self):
        """Test full pmap generation with GFE normalization"""
        from script.genpmap import gen_pmap
        import tempfile
        
        # Create minimal settings
        setting_general = {"name": "test"}
        setting_input = {
            "protein": {"pdb": str(self.ref_struct)},
            "probe": {"cid": "TEST"}
        }
        setting_pmap = {
            "snapshot": "1-10",
            "valid_dist": 18.0,
            "map_size": 80,
            "normalization": "GFE",
            "maps": [{"suffix": "test", "selector": "!@H*"}]
        }
        
        # Test pmap generation
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = gen_pmap(
                Path(tmpdir),
                setting_general,
                setting_input,
                setting_pmap,
                Path(tmpdir) / "traj.nc",
                Path(tmpdir) / "top.parm7"
            )
            self.assertTrue(len(paths) > 0)
            self.assertTrue(any("GFE_" in path for path in paths))

    def test_gen_pmap_others_normalization(self):
        """Test full pmap generation with total normalization"""
        from script.genpmap import gen_pmap
        import tempfile
        
        # Create minimal settings
        setting_general = {"name": "test"}
        setting_input = {
            "protein": {"pdb": str(self.ref_struct)},
            "probe": {"cid": "TEST"}
        }
        setting_pmap = {
            "snapshot": "1-10",
            "valid_dist": 18.0,
            "map_size": 80,
            "normalization": "total",
            "maps": [{"suffix": "test", "selector": "!@H*"}]
        }
        
        # Test pmap generation
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = gen_pmap(
                Path(tmpdir),
                setting_general,
                setting_input,
                setting_pmap,
                Path(tmpdir) / "traj.nc",
                Path(tmpdir) / "top.parm7"
            )
            self.assertTrue(len(paths) > 0)
