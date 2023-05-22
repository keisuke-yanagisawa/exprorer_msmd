from unittest import TestCase
from script.utilities.Bio import PDB


class TestPDB(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestPDB, self).__init__(*args, **kwargs)
        self.pdbfile = "script/utilities/Bio/test_data/PDB/7m67.pdb"

    def test_multi_model_pdb_reader(self):
        reader = PDB.MultiModelPDBReader(self.pdbfile)
        models = [model for model in reader]
        self.assertEqual(len(models), 10)

    def test_multi_model_pdb_reader_get_model_out_of_range(self):
        reader = PDB.MultiModelPDBReader(self.pdbfile)
        with self.assertRaises(IndexError):
            reader.get_model(10)

    def test_pdb_io_helper(self):
        reader = PDB.MultiModelPDBReader(self.pdbfile)
        models = [model for model in reader]

        outputfile = "script/utilities/Bio/test_data/PDB/output.pdb"
        writer = PDB.PDBIOhelper(outputfile)
        for model in models:
            writer.save(model)
        writer.close()  # This is important, otherwise the unittest will hang up

        reader = PDB.MultiModelPDBReader(outputfile)
        models = [model for model in reader]
        self.assertEqual(len(models), 10)

    def test_open_gziped_pdb_with_multi_model_pdb_reader(self):
        reader = PDB.MultiModelPDBReader("script/utilities/Bio/test_data/PDB/7m67.pdb.gz")
        models = [model for model in reader]
        self.assertEqual(len(models), 10)

    def test_get_structure_with_gziped_pdb(self):
        structure = PDB.get_structure("script/utilities/Bio/test_data/PDB/7m67.pdb.gz")
        self.assertEqual(len(structure), 10)
