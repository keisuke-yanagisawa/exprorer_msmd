from unittest import TestCase

import numpy as np

from script import alignresenv
from script.utilities.Bio import PDB as uPDB


class TestAlignResEnv(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestAlignResEnv, self).__init__(*args, **kwargs)
        self.single_atom_struct = uPDB.get_structure("script/test_data/singleatom.pdb")
        self.two_atoms_struct = uPDB.get_structure("script/test_data/twoatoms.pdb")
        self.two_models_struct = uPDB.get_structure("script/test_data/twomodels.pdb")
        self.no_atom_struct = uPDB.get_structure("script/test_data/noatom.pdb")

    def test_single_structure(self):
        alignresenv.align_res_env(self.single_atom_struct, self.single_atom_struct[0], "VAL")

    def test_single_structure_with_multiple_models(self):
        alignresenv.align_res_env(self.two_models_struct, self.two_models_struct[0], "PRO")

    def test_structure_with_no_atom_model(self):
        with self.assertRaises(ValueError):
            alignresenv.align_res_env(self.no_atom_struct, self.single_atom_struct[0], "VAL")

    def test_Reference_does_not_have_any_atom(self):
        with self.assertRaises(ValueError):
            alignresenv.align_res_env(self.two_models_struct, self.no_atom_struct[0], "UNK", focused=[" CA "])

    def test_focused_stable(self):
        alignresenv.align_res_env(
            self.two_models_struct, self.two_models_struct[0], "PRO", focused=[" CA ", " N  ", " C  "]
        )

    def test_focused_unstable(self):
        alignresenv.align_res_env(self.two_models_struct, self.two_models_struct[0], "PRO", focused=[" CA "])

    def test_focused_is_not_in_struct(self):
        with self.assertRaises(ValueError):
            alignresenv.align_res_env(self.two_models_struct, self.two_models_struct[0], "PRO", focused=[" XX "])

    def test_resn_not_in_structs(self):
        with self.assertRaises(ValueError):
            alignresenv.align_res_env(self.two_models_struct, self.two_models_struct[0], "UNK", focused=[" CA "])

    def test_aligned_coords(self):
        struct = alignresenv.align_res_env(self.two_atoms_struct, self.single_atom_struct[0], "VAL", focused=[" CA "])
        # TODO: placement of two_atoms_struct is the same from input to output
        #       because the position of CA atom is the same in both structures
        expected_coord = [[-0.603, 65.642, 77.183], [-0.883, 67.005, 76.630]]
        aligned_coord = [a.get_coord() for a in struct.get_atoms()]
        np.testing.assert_array_almost_equal(expected_coord, aligned_coord, decimal=3)
