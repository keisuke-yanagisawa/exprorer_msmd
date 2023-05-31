from unittest import TestCase
from script.resenv import resenv
import tempfile
from script.utilities.Bio import PDB as uPDB
import numpy as np


class TestResenv(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestResenv, self).__init__(*args, **kwargs)
        self.gridfile = "script/test_data/maxpmap_for_resenv.dx"
        self.trajectories = ["script/test_data/trajectory_for_resenv.pdb"]
        self.expected_resenv_pdb = "script/test_data/resenv_expected.pdb"
        self.resn = "A11"
        _, self.outputpdb = tempfile.mkstemp(suffix=".pdb")

    def test_normal_case(self):
        resenv(self.gridfile, self.trajectories, self.resn, self.outputpdb,
               threshold=0.001)

        struct = uPDB.get_structure(self.outputpdb)
        expected_struct = uPDB.get_structure(self.expected_resenv_pdb)
        self.assertEqual(len(struct), len(expected_struct))
        for model, expected_model in zip(struct, expected_struct):  # type: ignore
            np.testing.assert_array_almost_equal(uPDB.get_attr(model, "coord"), uPDB.get_attr(expected_model, "coord"))

    def test_grid_file_is_not_found(self):
        with self.assertRaises(FileNotFoundError):
            resenv("NOTHING", self.trajectories, self.resn, self.outputpdb)

    def test_ipdb_file_is_not_found(self):
        with self.assertRaises(FileNotFoundError):
            resenv(self.gridfile, ["NOTHING"], self.resn, self.outputpdb)

    def test_one_of_ipdb_file_is_not_found(self):
        with self.assertRaises(FileNotFoundError):
            resenv(self.gridfile, [self.trajectories[0], "NOTHING"], self.resn, self.outputpdb)

    def test_invalid_resn(self):
        with self.assertRaises(ValueError):
            resenv(self.gridfile, self.trajectories, "INVALID_RESN", self.outputpdb)

    def test_output_file_cannot_be_created(self):
        with self.assertRaises(FileNotFoundError):
            resenv(self.gridfile, self.trajectories, self.resn, "/INVALID/PATH/TO/OUTPUT")

    def test_no_output_structure(self):
        pass
        # with self.assertRaises(FileNotFoundError):
        #     resenv(self.gridfile, self.trajectories, "INVALID_RESN", None)
