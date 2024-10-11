import tempfile
from pathlib import Path
from unittest import TestCase

from script.utilities.Bio import PDB


class TestPDB(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestPDB, self).__init__(*args, **kwargs)
        self.pdbfile = Path("script/utilities/Bio/test_data/PDB/7m67.pdb")
        self.gziped_pdbfile = Path("script/utilities/Bio/test_data/PDB/7m67.pdb.gz")

    def test_multi_model_pdb_reader(self):
        reader = PDB.MultiModelPDBReader(str(self.pdbfile))
        models = [model for model in reader]
        self.assertEqual(len(models), 10)

    def test_multi_model_pdb_reader_get_model_out_of_range(self):
        reader = PDB.MultiModelPDBReader(str(self.pdbfile))
        with self.assertRaises(IndexError):
            reader.get_model(10)

    def test_multi_model_pdb_reader_not_file_found(self):
        with self.assertRaises(FileNotFoundError):
            PDB.MultiModelPDBReader("INVALID_PATH")
        pass

    def test_pdb_io_helper(self):
        reader = PDB.MultiModelPDBReader(str(self.pdbfile))
        models = [model for model in reader]
        tmp_output_pdb = Path(tempfile.mkstemp(suffix=".pdb")[1])
        writer = PDB.PDBIOhelper(tmp_output_pdb)
        for model in models:
            writer.save(model)
        writer.close()  # This is important, otherwise the unittest will hang up

        reader = PDB.MultiModelPDBReader(str(tmp_output_pdb))
        models = [model for model in reader]
        self.assertEqual(len(models), 10)

    def test_open_gziped_pdb_with_multi_model_pdb_reader(self):
        reader = PDB.MultiModelPDBReader(str(self.gziped_pdbfile))
        models = [model for model in reader]
        self.assertEqual(len(models), 10)

    def test_get_structure_with_gziped_pdb(self):
        structure = PDB.get_structure(self.gziped_pdbfile)
        self.assertEqual(len(structure), 10)
