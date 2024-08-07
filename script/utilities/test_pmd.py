import tempfile
from unittest import TestCase

from script.utilities import pmd


class TestConversion(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestConversion, self).__init__(*args, **kwargs)
        self.topfile = "script/utilities/test_data/pmd/system.parm7"
        self.xyzfile = "script/utilities/test_data/pmd/system.rst7"
        self.expected_topfile = "script/utilities/test_data/pmd/system.top"
        self.expected_grofile = "script/utilities/test_data/pmd/system.gro"

    def test_run_convert(self):
        _, out_top = tempfile.mkstemp(suffix=".top")
        _, out_gro = tempfile.mkstemp(suffix=".gro")
        pmd.convert(self.topfile, out_top, self.xyzfile, out_gro)

        # Ignore the first 15 lines including the timestamp and command to run test
        self.assertEqual(open(out_top).readlines()[15:], open(self.expected_topfile).readlines()[15:])
        self.assertEqual(open(out_gro).read(), open(self.expected_grofile).read())

    def test_run_convert_top_only(self):
        _, out_top = tempfile.mkstemp(suffix=".top")
        pmd.convert(self.topfile, out_top)

        # Ignore the first 15 lines including the timestamp and command to run test
        self.assertEqual(open(out_top).readlines()[15:], open(self.expected_topfile).readlines()[15:])

    def __del__(self):
        pass
        # os.system("rm -rf script/utilities/test_data/pmd/output")
