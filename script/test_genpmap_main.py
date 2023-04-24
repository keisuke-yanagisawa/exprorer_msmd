from unittest import TestCase
from genpmap_main import mask_generator
from utilities import util
import gridData
from utilities.Bio import PDB as uPDB
from utilities import GridUtil
import numpy as np


class TestMaskGenerator(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestMaskGenerator, self).__init__(*args, **kwargs)
        self.grid = gridData.Grid()
        self.grid.load("test_data/small_grid.dx")
        self.ref_struct = "test_data/tripeptide.pdb"

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
