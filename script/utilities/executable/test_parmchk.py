import tempfile
from pathlib import Path
from unittest import TestCase

from script.utilities.executable.parmchk import Parmchk


class TestParmchk(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestParmchk, self).__init__(*args, **kwargs)
        __testdata_dir = Path("script/utilities/executable/test_data")
        self.probe_mol2 = __testdata_dir / "A11.mol2"
        self.expected_probe_frcmod = __testdata_dir / "A11.frcmod"

    def test_parmchk(self):
        tmp_frcmod = Path(tempfile.mkstemp(suffix=".frcmod")[1])
        parmchk = Parmchk()
        parmchk.set(mol2=self.probe_mol2, at="gaff2")
        parmchk.run(frcmod=tmp_frcmod)

        expected = open(self.expected_probe_frcmod, "r").read()
        actual = open(tmp_frcmod, "r").read()
        self.assertEqual(expected, actual)

    def test_file_does_not_exist(self):
        parmchk = Parmchk()
        with self.assertRaises(FileNotFoundError):
            parmchk.set(mol2=Path("NOTHING.mol2"), at="gaff2")
            parmchk.run()

    def test_at_is_not_supported(self):
        parmchk = Parmchk()
        with self.assertRaises(ValueError):
            parmchk.set(mol2=self.probe_mol2, at="NOT_SUPPORTED_ID")  # type: ignore
            parmchk.run()
