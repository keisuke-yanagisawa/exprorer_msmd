from unittest import TestCase
from generate_msmd_system import calculate_boxsize, generate_msmd_system
from utilities import util
import os


class TestBoxSizeCalculation(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestBoxSizeCalculation, self).__init__(*args, **kwargs)

    def test_calculate_boxsize(self):
        rstfile = "test_data/tripeptide.rst7"
        box_size = calculate_boxsize(rstfile)
        expected = 16.4927710
        self.assertAlmostEqual(box_size, expected)

    def test_calculate_boxsize_error(self):
        pdbfile = "test_data/tripeptide.pdb"
        with self.assertRaises(ValueError):
            calculate_boxsize(pdbfile)


class TestGenerateMsmdSystem(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestGenerateMsmdSystem, self).__init__(*args, **kwargs)

    def test_generate_msmd_system(self):
        setting = util.parse_yaml("test_data/setting.yaml")
        parm7, rst7 = generate_msmd_system(setting, seed=1)

        # first line contains the date and time, so we skip it for comparison
        parm7_str = open(parm7).readlines()[1:]
        expected_parm7 = open("test_data/tripeptide_A11.parm7").readlines()[1:]
        self.assertEqual(parm7_str, expected_parm7)

        rst7_str = open(rst7).read()
        expected_rst7 = open("test_data/tripeptide_A11.rst7").read()
        self.assertEqual(rst7_str, expected_rst7)
