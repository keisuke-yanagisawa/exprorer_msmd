import tempfile
from pathlib import Path
from unittest import TestCase
from unittest.mock import MagicMock, patch

from script.utilities.Bio import PDB as uPDB
from script.utilities.executable.cpptraj import Cpptraj

# TODO: Add tests for cpptraj


class TestCpptraj(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestCpptraj, self).__init__(*args, **kwargs)
        self.cpptraj = Cpptraj()
        self.trajectory_path = Path("script/utilities/executable/test_data/cpptraj/trajectory.xtc")
        self.topology_path = Path("script/utilities/executable/test_data/cpptraj/topology.top")
        self.ref_struct_path = Path("script/utilities/executable/test_data/cpptraj/inputprotein.pdb")
        self.probe_id = "A11"
        self.maps = [{"suffix": "nVH", "selector": "(!@VIS)&(!@H*)"}]
        self.box_size = 80
        self.box_center = uPDB.get_structure(self.ref_struct_path).center_of_mass()

    def test_cpptraj_execution(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cpptraj_obj = Cpptraj()
            cpptraj_obj.set(self.topology_path, self.trajectory_path, self.ref_struct_path, self.probe_id)
            cpptraj_obj.run(
                basedir=Path(tmpdir),
                prefix="TEST",
                box_size=self.box_size,
                box_center=self.box_center,
                traj_start=1,
                traj_stop=5,
                traj_offset=1,
                maps=self.maps,
            )

    def test_is_cpptraj_output_indentical(self):
        pass

    @patch("script.utilities.executable.cpptraj.Command")
    def test_cpptraj_error_propagates(self, mock_command_cls):
        """Characterization test: when the cpptraj subprocess fails,
        the exception is re-raised to the caller (not swallowed).

        Pinning current behavior before refactoring the exception handler.
        """
        # Make the first Command (cpptraj invocation) fail; a second Command
        # may be constructed for the "cat" debug call but is not expected to
        # be .run()
        failing_command = MagicMock()
        failing_command.run.side_effect = RuntimeError("cpptraj failed")
        # Subsequent Command(...) calls (e.g. the `cat` debug one) are ok
        mock_command_cls.side_effect = [failing_command, MagicMock()]

        with tempfile.TemporaryDirectory() as tmpdir:
            cpptraj_obj = Cpptraj()
            cpptraj_obj.set(self.topology_path, self.trajectory_path, self.ref_struct_path, self.probe_id)
            with self.assertRaises(RuntimeError, msg="cpptraj failure should propagate"):
                cpptraj_obj.run(
                    basedir=Path(tmpdir),
                    prefix="TEST",
                    box_size=self.box_size,
                    box_center=self.box_center,
                    traj_start=1,
                    traj_stop=5,
                    traj_offset=1,
                    maps=self.maps,
                )
        # The main cpptraj command must have been attempted
        failing_command.run.assert_called_once()
