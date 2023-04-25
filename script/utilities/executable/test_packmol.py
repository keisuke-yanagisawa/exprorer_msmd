from unittest import TestCase
from script.utilities.executable.packmol import Packmol


class TestPackmol(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestPackmol, self).__init__(*args, **kwargs)

    def test_run_packmol(self):
        packmol = Packmol()
        packmol.set(protein_pdb="script/utilities/executable/test_data/tripeptide.pdb",
                    cosolv_pdb="script/utilities/executable/test_data/A11.pdb",
                    box_size=20,
                    molar=0.1)
        packmol.run()

    def test_zero_molar(self):
        with self.assertWarns(RuntimeWarning):
            packmol = Packmol()
            packmol.set(protein_pdb="script/utilities/executable/test_data/tripeptide.pdb",
                        cosolv_pdb="script/utilities/executable/test_data/A11.pdb",
                        box_size=20,
                        molar=0)
            packmol.run()

    def test_extremely_low_molar(self):
        packmol = Packmol()
        packmol.set(protein_pdb="script/utilities/executable/test_data/tripeptide.pdb",
                    cosolv_pdb="script/utilities/executable/test_data/A11.pdb",
                    box_size=20,
                    molar=1e-10)
        packmol.run()

    def test_high_molar(self):
        with self.assertWarns(RuntimeWarning):
            packmol = Packmol()
            packmol.set(protein_pdb="script/utilities/executable/test_data/tripeptide.pdb",
                        cosolv_pdb="script/utilities/executable/test_data/A11.pdb",
                        box_size=20,
                        molar=10)
            packmol.run()

    def test_extremely_high_molar(self):
        with self.assertRaises(RuntimeError):
            packmol = Packmol()
            packmol.set(protein_pdb="script/utilities/executable/test_data/tripeptide.pdb",
                        cosolv_pdb="script/utilities/executable/test_data/A11.pdb",
                        box_size=20,
                        molar=1e10)
            packmol.run()

    def test_no_protein(self):
        with self.assertRaises(FileNotFoundError):
            packmol = Packmol()
            packmol.set(protein_pdb="NOTHING.pdb",
                        cosolv_pdb="script/utilities/executable/test_data/A11.pdb",
                        box_size=20,
                        molar=0.1)
            packmol.run()

    def test_no_cosolvent(self):
        with self.assertRaises(FileNotFoundError):
            packmol = Packmol()
            packmol.set(protein_pdb="script/utilities/executable/test_data/tripeptide.pdb",
                        cosolv_pdb="NOTHING.pdb",
                        box_size=20,
                        molar=0.1)
            packmol.run()

    def test_cosolvent_is_not_pdb(self):
        with self.assertRaises(RuntimeError):
            packmol = Packmol()
            packmol.set(protein_pdb="script/utilities/executable/test_data/tripeptide.pdb",
                        cosolv_pdb="script/utilities/executable/test_data/A11.mol2",
                        box_size=20,
                        molar=0.1)
            packmol.run()

    def test_cosolvent_file_contains_non_cosolvent(self):
        with self.assertRaises(RuntimeError):
            packmol = Packmol()
            packmol.set(protein_pdb="script/utilities/executable/test_data/tripeptide.pdb",
                        cosolv_pdb="script/utilities/executable/test_data/tripeptide.pdb",
                        box_size=20,
                        molar=0.1)
            packmol.run()

    def test_protein_file_contains_non_protein(self):
        with self.assertRaises(RuntimeError):
            packmol = Packmol()
            packmol.set(protein_pdb="script/utilities/executable/test_data/A11.pdb",
                        cosolv_pdb="script/utilities/executable/test_data/A11.pdb",
                        box_size=20,
                        molar=0.1)
            packmol.run()
