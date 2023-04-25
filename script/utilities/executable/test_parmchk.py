from unittest import TestCase
from script.utilities.executable.parmchk import Parmchk


class TestParmchk(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestParmchk, self).__init__(*args, **kwargs)

    def test_parmchk(self):
        parmchk = Parmchk()
        parmchk.set(mol2="script/utilities/executable/test_data/A11.mol2",
                    at="gaff2")
        parmchk.run(frcmod="script/utilities/executable/test_data/output/A11.frcmod")
        with open("script/utilities/executable/test_data/A11.frcmod", "r") as f:
            expected = f.read()
        with open("script/utilities/executable/test_data/output/A11.frcmod", "r") as f:
            actual = f.read()
        self.assertEqual(expected, actual)

    def test_file_does_not_exist(self):
        parmchk = Parmchk()
        with self.assertRaises(FileNotFoundError):
            parmchk.set(mol2="NOTHING.mol2",
                        at="gaff2")
            parmchk.run()

    def test_at_is_not_supported(self):
        parmchk = Parmchk()
        with self.assertRaises(ValueError):
            parmchk.set(mol2="script/utilities/executable/test_data/A11.mol2",
                        at="NOT_SUPPORTED_ID")
            parmchk.run()
