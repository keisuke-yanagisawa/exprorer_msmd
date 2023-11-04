from unittest import TestCase
from script.utilities.executable.cpptraj import Cpptraj
from script.utilities.Bio import PDB as uPDB
# TODO: Add tests for cpptraj


class TestCpptraj(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestCpptraj, self).__init__(*args, **kwargs)
        self.cpptraj = Cpptraj()
        self.trajectory_path = "script/utilities/executable/test_data/cpptraj/trajectory.xtc"
        self.topology_path = "script/utilities/executable/test_data/cpptraj/topology.top"
        self.ref_struct_path = "script/utilities/executable/test_data/cpptraj/inputprotein.pdb"
        self.probe_id = "A11"
        self.maps = [{"suffix": "nVH", "selector": "(!@VIS)&(!@H*)"}]
        self.box_size = 80
        self.box_center = uPDB.get_structure(self.ref_struct_path).center_of_mass()

    def test_cpptraj_execution(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            cpptraj_obj = Cpptraj()
            cpptraj_obj.set(self.topology_path, self.trajectory_path, self.ref_struct_path, self.probe_id)
            cpptraj_obj.run(
                basedir=tmpdir,
                prefix="TEST",
                box_size=self.box_size,
                box_center=self.box_center,
                traj_start=1,
                traj_stop=5,
                traj_offset=1,
                maps=self.maps
            )

    def test_is_cpptraj_output_indentical(self):
        pass
