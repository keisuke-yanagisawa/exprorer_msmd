from pathlib import Path
from unittest import TestCase

import gridData
import numpy as np

from script.resenv import resenv
from script.utilities.Bio import PDB as uPDB


class TestResenv(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestResenv, self).__init__(*args, **kwargs)
        __testdata_dir = Path("script/test_data")
        self.gridfile = __testdata_dir / "maxpmap_for_resenv.dx"
        self.grid = gridData.Grid(self.gridfile)
        self.trajectory_file = __testdata_dir / "trajectory_for_resenv.pdb"
        self.trajectory = uPDB.MultiModelPDBReader(str(self.trajectory_file))
        self.expected_resenv_pdb = __testdata_dir / "resenv_expected.pdb"
        self.resn = "A11"

    def test_normal_case(self):
        struct = resenv(self.grid, self.trajectory, self.resn, [" CB "], threshold=0.001)
        expected_struct = uPDB.get_structure(self.expected_resenv_pdb)
        self.assertEqual(len(struct), len(expected_struct))
        for model, expected_model in zip(struct, expected_struct):  # type: ignore
            np.testing.assert_array_almost_equal(uPDB.get_attr(model, "coord"), uPDB.get_attr(expected_model, "coord"))

    def test_invalid_resn(self):
        with self.assertRaises(ValueError):
            resenv(self.grid, self.trajectory, "INVALID_RESN", [" CB "])

    def test_no_output_structure(self):
        pass
        # with self.assertRaises(FileNotFoundError):
        #     resenv(self.grid, self.trajectories, "INVALID_RESN", None)
