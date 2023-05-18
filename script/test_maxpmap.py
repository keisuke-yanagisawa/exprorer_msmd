from unittest import TestCase
from script.maxpmap import grid_max
import gridData

class TestMaskGenerator(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestMaskGenerator, self).__init__(*args, **kwargs)
        self.grid1 = gridData.Grid("script/test_data/pmap1.dx")
        self.grid2 = gridData.Grid("script/test_data/pmap2.dx")

    def test_maxpmap(self):
        maxpmap = grid_max([self.grid1, self.grid2])
        self.assertAlmostEqual(self.grid1.grid[1,5,2], maxpmap.grid[1,5,2]) # -1
        self.assertAlmostEqual(self.grid1.grid[38,59,27], maxpmap.grid[38,59,27]) # 0
        self.assertAlmostEqual(self.grid1.grid[44,20,32], maxpmap.grid[44,20,32]) # not -1 & not 0 for grid1, grid2[44,20,32] = 0
