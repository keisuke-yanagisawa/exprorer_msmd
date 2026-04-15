from pathlib import Path
from unittest import TestCase

import numpy as np

from script import alignresenv
from script.utilities.Bio import PDB as uPDB


class TestAlignResEnv(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestAlignResEnv, self).__init__(*args, **kwargs)
        __testdata_dir = Path("script/test_data")

        self.single_atom_struct = uPDB.get_structure(__testdata_dir / "singleatom.pdb")
        self.two_atoms_struct = uPDB.get_structure(__testdata_dir / "twoatoms.pdb")
        self.two_models_struct = uPDB.get_structure(__testdata_dir / "twomodels.pdb")
        self.no_atom_struct = uPDB.get_structure(__testdata_dir / "noatom.pdb")

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

    def test_two_models_stable_focus_coords(self):
        """Characterization test: verify coords of all models after stable (3-atom) alignment.

        Snapshot captured before refactoring of alignresenv tempfile handling
        to ensure the refactored implementation produces identical output.
        """
        struct = alignresenv.align_res_env(
            self.two_models_struct,
            self.two_models_struct[0],
            "PRO",
            focused=[" CA ", " N  ", " C  "],
        )
        # Expected coords per model: model0 and model1 both have (N, CA)
        expected = [
            [4.524, 9.887, -0.667],
            [5.918, 10.123, -0.175],
            [4.524, 9.887, -0.667],
            [5.918, 10.123, -0.175],
        ]
        actual = [a.get_coord() for a in struct.get_atoms()]
        np.testing.assert_array_almost_equal(expected, actual, decimal=3)

    def test_two_models_unstable_focus_coords(self):
        """Characterization test: verify coords after unstable (single-atom) alignment.

        With only CA atom as focus, non-focused atoms (N) may differ between
        models due to rotational ambiguity. Snapshot captured before refactoring.
        """
        struct = alignresenv.align_res_env(
            self.two_models_struct,
            self.two_models_struct[0],
            "PRO",
            focused=[" CA "],
        )
        expected = [
            [4.524, 9.887, -0.667],
            [5.918, 10.123, -0.175],
            [5.259, 11.002, -1.193],
            [5.918, 10.123, -0.175],
        ]
        actual = [a.get_coord() for a in struct.get_atoms()]
        np.testing.assert_array_almost_equal(expected, actual, decimal=3)
