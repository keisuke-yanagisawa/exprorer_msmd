from unittest import TestCase
from script.resenv import resenv
import tempfile


class TestResenvMain(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestResenvMain, self).__init__(*args, **kwargs)
        self.gridfile = "script/test_data/maxpmap_for_resenv.dx"
        self.trajectories = ["script/test_data/trajectory_for_resenv.pdb"]
        self.expected_resenv_pdb = "script/test_data/resenv_expected.pdb"
        self.resn = "A11"

    def test_resenv(self):
        _, tmppdb = tempfile.mkstemp(suffix=".pdb")
        resenv(self.gridfile, self.trajectories, self.resn, tmppdb,
               threshold=0.001)

        resenv_str = open(tmppdb).read()
        expected = open(self.expected_resenv_pdb).read()
        self.assertEqual(resenv_str, expected)
