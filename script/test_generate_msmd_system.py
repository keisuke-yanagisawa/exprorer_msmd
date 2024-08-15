from pathlib import Path
from typing import Literal
from unittest import TestCase

from script.generate_msmd_system import _create_frcmod, calculate_boxsize, generate_msmd_system
from script.utilities import util


class TestBoxSizeCalculation(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestBoxSizeCalculation, self).__init__(*args, **kwargs)
        __testdata_dir = Path("script/test_data")
        self.rstfile = __testdata_dir / "tripeptide.rst7"
        self.pdbfile = __testdata_dir / "tripeptide.pdb"

    def test_calculate_boxsize(self):
        box_size = calculate_boxsize(self.rstfile)
        expected = 16.4927710
        self.assertAlmostEqual(box_size, expected)

    def test_calculate_boxsize_error(self):
        with self.assertRaises(ValueError):
            calculate_boxsize(self.pdbfile)


class TestGenerateMsmdSystem(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestGenerateMsmdSystem, self).__init__(*args, **kwargs)
        __testdata_dir = Path("script/test_data")
        self.setting_file = __testdata_dir / "setting.yaml"
        self.expected_parm7_file = __testdata_dir / "tripeptide_A11.parm7"
        self.expected_rst7_file = __testdata_dir / "tripeptide_A11.rst7"

    def test_generate_msmd_system(self):
        settings = util.parse_yaml(self.setting_file)
        parm7, rst7 = generate_msmd_system(settings, seed=1)

        # first line contains the date and time, so we skip it for comparison
        parm7_str = open(parm7).readlines()[1:]
        expected_parm7 = open(self.expected_parm7_file).readlines()[1:]
        self.assertEqual(parm7_str, expected_parm7)

        rst7_str = open(rst7).read()
        expected_rst7 = open(self.expected_rst7_file).read()
        self.assertEqual(rst7_str, expected_rst7)


class TestCreateFrcmod(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestCreateFrcmod, self).__init__(*args, **kwargs)
        __testdata_dir = Path("script/test_data")
        self.mol2file = __testdata_dir / "A11.mol2"
        self.pdbfile = __testdata_dir / "A11.pdb"
        self.expected_frcmod = __testdata_dir / "A11.frcmod"
        self.atomtype: Literal["gaff", "gaff2"] = "gaff2"
        self.invalid_path = Path("INVALID")

    def test_create_frcmod(self):
        frcmod_path = _create_frcmod(self.mol2file, self.atomtype)

        frcmod_str = open(frcmod_path).read()
        expected_frcmod = open(self.expected_frcmod).read()
        self.assertEqual(frcmod_str, expected_frcmod)

    def test_invalid_atomtype(self):
        with self.assertRaises(ValueError):
            _create_frcmod(self.mol2file, self.invalid_path)  # type: ignore

    def test_mol2file_does_not_exist(self):
        with self.assertRaises(FileNotFoundError):
            _create_frcmod(self.invalid_path, self.atomtype)  # type: ignore

    def test_invalid_mol2file(self):
        with self.assertRaises(ValueError):
            _create_frcmod(self.pdbfile, self.atomtype)
