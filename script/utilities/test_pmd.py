import os
from unittest import TestCase
from script.utilities import pmd


class TestConversion(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestConversion, self).__init__(*args, **kwargs)
        self.topfile = "script/utilities/test_data/pmd/system.parm7"
        self.xyzfile = "script/utilities/test_data/pmd/system.rst7"

    def test_run_convert(self):
        os.system("mkdir -p script/utilities/test_data/pmd/output")
        out_top = "script/utilities/test_data/pmd/output/system.top"
        out_gro = "script/utilities/test_data/pmd/output/system.gro"
        pmd.convert(self.topfile, out_top, self.xyzfile, out_gro)

        expected_top = "script/utilities/test_data/pmd/system.top"
        expected_gro = "script/utilities/test_data/pmd/system.gro"

        self.assertEqual(open(out_top).readlines()[6:],
                         open(expected_top).readlines()[6:])  # Ignore the first 6 lines including the timestamp
        self.assertEqual(open(out_gro).read(),
                         open(expected_gro).read())

    def test_run_convert_top_only(self):
        os.system("mkdir -p script/utilities/test_data/pmd/output")
        out_top = "script/utilities/test_data/pmd/output/system.top"
        pmd.convert(self.topfile, out_top)

        expected_top = "script/utilities/test_data/pmd/system.top"

        self.assertEqual(open(out_top).readlines()[6:],
                         open(expected_top).readlines()[6:])  # Ignore the first 6 lines including the timestamp

    def __del__(self):
        pass
        # os.system("rm -rf script/utilities/test_data/pmd/output")
