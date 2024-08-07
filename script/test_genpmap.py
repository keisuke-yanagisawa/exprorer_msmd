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
        pass

    def test_convert_to_proba_GFE_normalization(self):
        pass

    def test_convert_to_proba_others_normalization(self):
        pass

    def test_convert_to_gfe(self):
        pass

    def test_convert_to_pmap_snapshot_normalization(self):
        pass

    def test_convert_to_pmap_others_normalization(self):
        pass

    def test_parse_snapshot_setting(self):
        pass

    def test_gen_pmap_snapshot_normalization(self):
        pass

    def test_gen_pmap_gfe_normalization(self):
        pass

    def test_gen_pmap_others_normalization(self):
        pass
