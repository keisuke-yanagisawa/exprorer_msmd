from unittest import TestCase
from script.utilities.GridUtil import gen_distance_grid
import gridData


class TestGenDistanceGrid(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestGenDistanceGrid, self).__init__(*args, **kwargs)
        self.grid = gridData.Grid()
        self.grid.load("script/test_data/small_grid.dx")
        self.input_pdb = "script/test_data/tripeptide.pdb"

    def test_gen_distance_grid(self):
        g = gen_distance_grid(self.grid, self.input_pdb)
        self.assertEqual(g.grid.shape, (2, 2, 2))
        self.assertAlmostEqual(g.grid[0, 0, 0], 17.72341, places=3)
