import tempfile
import os
from subprocess import getoutput as gop
from scipy import constants
from utilities import const
from utilities.executable.execute import Command
import shutil
import jinja2
from utilities.pmd import convert as pmd_convert
from ..logger import logger

class Cpptraj(object):
    def __init__(self, debug=False):
        self.exe = os.getenv("CPPTRAJ", "cpptraj")
        self.debug = debug

    def _gen_parm7(self):
        _, self.parm7 = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PARM7)
        self.parm7, _ = pmd_convert(self.topology, self.parm7)

    def set(self, topology, trajectory, ref_struct, probe_id):
        self.topology = topology
        self.trajectory = trajectory
        self.ref_struct = ref_struct
        self.probe_id = probe_id
        return self

    def run(self, basedir, prefix, box_size=80, interval=1, 
            traj_start=1, traj_stop="last", traj_offset=1):
        self.basedir = basedir
        self.prefix = prefix
        self.voxel = [box_size, interval] * 3 # x, y, z
        self.frame_info = [traj_start, traj_stop, traj_offset]

        self._gen_parm7()

        _, self.inp = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)
        data = {
          "basedir": self.basedir,
          "top" : self.parm7,
          "traj": self.trajectory,
          "cid": self.probe_id,
          "frame_info": " ".join([str(n) for n in self.frame_info]),
          "ref": self.ref_struct,
          "map_voxel": " ".join([str(n) for n in self.voxel]),
          "prefix": self.prefix
        }

        env = jinja2.Environment(loader=jinja2.FileSystemLoader(f"{os.path.dirname(__file__)}/template"))
        template = env.get_template("cpptraj_pmap.in")
        with open(self.inp, "w") as fout:
            fout.write(template.render(data))
            print(template.render(data))
        command = Command(f"{self.exe} < {self.inp}")
        if self.debug:
          print(command)
        print(command.run())

        self.grids = [
            f"{self.basedir}/{self.prefix}_nV.dx",
            f"{self.basedir}/{self.prefix}_nVH.dx",
            f"{self.basedir}/{self.prefix}_V.dx",
            f"{self.basedir}/{self.prefix}_O.dx",
        ]

        return self

    def __del__(self):
      logger.debug(f"Cpptraj.inp: {self.inp}")
      logger.debug(f"Cpptraj.parm7: {self.parm7}")
      if hasattr(self, "inp"):
        if not self.debug:
          os.remove(self.inp)
      if hasattr(self, "parm7"):
        if not self.debug:
          os.remove(self.parm7)
