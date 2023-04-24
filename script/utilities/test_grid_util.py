from unittest import TestCase
from utilities.util import parse_yaml, expand_index
from utilities.GridUtil import gen_distance_grid
import gridData
from utilities.Bio import PDB as uPDB


class TestGenDistanceGrid(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestGenDistanceGrid, self).__init__(*args, **kwargs)
        self.grid = gridData.Grid()
        self.grid.load("test_data/small_grid.dx")

    def test_gen_distance_grid(self):
        g = gen_distance_grid(self.grid, "test_data/tripeptide.pdb")
        self.assertEqual(g.grid.shape, (2, 2, 2))
        self.assertAlmostEqual(g.grid[0, 0, 0], 17.72341, places=3)
