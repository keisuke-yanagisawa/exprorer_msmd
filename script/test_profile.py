from pathlib import Path
from unittest import TestCase

import numpy as np

from script import profile
from script.utilities.Bio import PDB as uPDB


class TestCreateResidueInteractionProfile(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestCreateResidueInteractionProfile, self).__init__(*args, **kwargs)
        __test_data = Path("script/test_data/")
        self.single_atom_struct = uPDB.get_structure(__test_data / "singleatom.pdb")
        self.two_atoms_struct = uPDB.get_structure(__test_data / "twoatoms.pdb")
        self.two_models_struct = uPDB.get_structure(__test_data / "twomodels.pdb")
        self.no_atom_struct = uPDB.get_structure(__test_data / "noatom.pdb")

    def test_normal_case(self):
        target_residue_atoms = [("VAL", " CA "), ("VAL", " N  ")]
        grid = profile.create_residue_interaction_profile(self.two_atoms_struct, target_residue_atoms)
        np.testing.assert_array_equal(grid.grid, [[[0, 1], [0, 0], [1, 0]]])
        np.testing.assert_array_almost_equal(grid.origin, [-0.5, 65.5, 76.5])

    def test_one_atom(self):

        grid = profile.create_residue_interaction_profile(self.two_atoms_struct, [("VAL", " CA ")])
        np.testing.assert_array_equal(grid.grid, [[[1]]])
        np.testing.assert_array_almost_equal(grid.origin, [-0.5, 67.5, 76.5])

    def test_multi_models(self):
        grid = profile.create_residue_interaction_profile(self.two_models_struct, [("PRO", " N  ")])
        np.testing.assert_array_equal(grid.grid, [[[1], [0]], [[0], [0]], [[0], [0]], [[0], [0]], [[0], [1]]])
        np.testing.assert_array_almost_equal(grid.origin, [4.5, 9.5, -0.5])

    def test_empty_model(self):
        self.assertEqual(len([m for m in self.no_atom_struct.get_models()]), 1)
        self.assertEqual(len([a for a in self.no_atom_struct.get_atoms()]), 0)
        with self.assertRaises(ValueError):
            profile.create_residue_interaction_profile(self.no_atom_struct, [("VAL", " CA ")])

    def test_no_residue(self):
        with self.assertRaises(ValueError):
            profile.create_residue_interaction_profile(self.two_atoms_struct, [("ALA", " CA ")])

    def test_no_atomname(self):
        with self.assertRaises(ValueError):
            profile.create_residue_interaction_profile(self.two_atoms_struct, [("VAL", " CB ")])
