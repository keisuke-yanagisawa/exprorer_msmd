import tempfile
from pathlib import Path
from unittest import TestCase

from script.utilities import pmd


class TestConversion(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestConversion, self).__init__(*args, **kwargs)
        __testdata_dir = Path("script/utilities/test_data")
        self.topfile = __testdata_dir / "pmd/system.parm7"
        self.xyzfile = __testdata_dir / "pmd/system.rst7"
        self.expected_topfile = __testdata_dir / "pmd/system.top"
        self.expected_grofile = __testdata_dir / "pmd/system.gro"

    def test_run_convert(self):
        out_top = Path(tempfile.mkstemp(suffix=".top")[1])
        out_gro = Path(tempfile.mkstemp(suffix=".gro")[1])
        pmd.convert(self.topfile, out_top, self.xyzfile, out_gro)

        # Ignore the first 15 lines including the timestamp and command to run test
        self.assertEqual(open(out_top).readlines()[15:], open(self.expected_topfile).readlines()[15:])
        self.assertEqual(open(out_gro).read(), open(self.expected_grofile).read())

    def test_run_convert_top_only(self):
        out_top = Path(tempfile.mkstemp(suffix=".top")[1])
        pmd.convert(self.topfile, out_top)

        # Ignore the first 15 lines including the timestamp and command to run test
        self.assertEqual(open(out_top).readlines()[15:], open(self.expected_topfile).readlines()[15:])

    def __del__(self):
        pass
        # os.system("rm -rf script/utilities/test_data/pmd/output")
