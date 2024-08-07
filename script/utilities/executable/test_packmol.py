import warnings
from pathlib import Path
from unittest import TestCase

from script.utilities.executable.packmol import Packmol


class TestPackmol(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestPackmol, self).__init__(*args, **kwargs)
        __testdata_dir = Path("script/utilities/executable/test_data")
        self.protein_pdb = __testdata_dir / "tripeptide.pdb"
        self.cosolv_pdb = __testdata_dir / "A11.pdb"
        self.cosolv_mol2 = __testdata_dir / "A11.mol2"

    def test_run_packmol(self):
        packmol = Packmol()
        packmol.set(protein_pdb=self.protein_pdb, cosolv_pdb=self.cosolv_pdb, box_size=20, molar=0.5)
        packmol.run(seed=1)

    def test_zero_molar(self):
        with self.assertWarns(RuntimeWarning):
            packmol = Packmol()
            packmol.set(protein_pdb=self.protein_pdb, cosolv_pdb=self.cosolv_pdb, box_size=20, molar=0)
            packmol.run(seed=1)

    def test_extremely_low_molar(self):
        packmol = Packmol()
        with self.assertWarns(RuntimeWarning):
            # the concentration is too low to add even single probe molecule
            packmol.set(protein_pdb=self.protein_pdb, cosolv_pdb=self.cosolv_pdb, box_size=20, molar=1e-10)
            packmol.run(seed=1)

    def test_high_molar(self):
        with self.assertWarns(RuntimeWarning):
            packmol = Packmol()
            packmol.set(protein_pdb=self.protein_pdb, cosolv_pdb=self.cosolv_pdb, box_size=20, molar=10)
            packmol.run(seed=1)

    def test_extremely_high_molar(self):
        with self.assertRaises(RuntimeError):
            packmol = Packmol()
            with warnings.catch_warnings():
                # this test case will raise a RuntimeWarning, but we don't care
                warnings.simplefilter("ignore", RuntimeWarning)
                packmol.set(protein_pdb=self.protein_pdb, cosolv_pdb=self.cosolv_pdb, box_size=20, molar=1e10)
            packmol.run(seed=1)

    def test_no_protein(self):
        with self.assertRaises(FileNotFoundError):
            packmol = Packmol()
            packmol.set(protein_pdb=Path("NOTHING.pdb"), cosolv_pdb=self.cosolv_pdb, box_size=20, molar=0.1)
            packmol.run(seed=1)

    def test_no_cosolvent(self):
        with self.assertRaises(FileNotFoundError):
            packmol = Packmol()
            packmol.set(protein_pdb=self.protein_pdb, cosolv_pdb=Path("NOTHING.pdb"), box_size=20, molar=0.1)
            packmol.run(seed=1)

    def test_cosolvent_is_not_pdb(self):
        with self.assertRaises(ValueError):
            packmol = Packmol()
            packmol.set(protein_pdb=self.protein_pdb, cosolv_pdb=self.cosolv_mol2, box_size=20, molar=0.1)
            packmol.run(seed=1)

    def test_cosolvent_file_contains_non_cosolvent(self):
        with self.assertRaises(RuntimeError):
            packmol = Packmol()
            packmol.set(protein_pdb=self.protein_pdb, cosolv_pdb=self.protein_pdb, box_size=20, molar=0.1)
            packmol.run(seed=1)

    def test_protein_file_contains_non_protein(self):
        with self.assertRaises(RuntimeError):
            packmol = Packmol()
            packmol.set(protein_pdb=self.cosolv_pdb, cosolv_pdb=self.cosolv_pdb, box_size=20, molar=0.1)
            packmol.run(seed=1)
