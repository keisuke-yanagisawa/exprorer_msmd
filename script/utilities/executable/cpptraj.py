import copy
import os
import tempfile
from pathlib import Path
from typing import Union

import jinja2
import numpy as np
import numpy.typing as npt

from .. import const
from ..logger import logger
from ..pmd import convert as pmd_convert
from .execute import Command


class Cpptraj(object):
    def __init__(self, debug: bool = False):
        self.exe: Path = Path(os.getenv("CPPTRAJ", "cpptraj"))
        self.debug = debug

    def _gen_parm7(self) -> None:
        self.parm7 = Path(tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PARM7)[1])
        self.parm7, _ = pmd_convert(self.topology, self.parm7)

    def set(self, topology: Path, trajectory: Path, ref_struct: Path, probe_id: str) -> "Cpptraj":
        self.topology = topology
        self.trajectory = trajectory
        self.ref_struct = ref_struct
        self.probe_id = probe_id

        return self

    def run(
        self,
        basedir: Path,
        prefix: str,
        box_center: npt.NDArray[np.float_] = np.array([0.0, 0.0, 0.0]),
        box_size: int = 80,
        interval: float = 1.0,
        traj_start: Union[str, int] = 1,
        traj_stop: Union[str, int] = "last",
        traj_offset: Union[str, int] = 1,
        maps: list = [{"suffix": "nVH", "selector": "(!@VIS)&(!@H*)"}],
    ):
        # TODO: input "maps" variable should be a read-only list (shared between threads)
        maps = copy.deepcopy(maps)
        self.basedir = basedir
        self.prefix = prefix
        self.voxel: list[Union[int, float]] = [box_size, interval] * 3  # x, y, z
        self.frame_info: tuple[Union[str, int], Union[str, int], Union[str, int]] = (traj_start, traj_stop, traj_offset)
        for i in range(len(maps)):
            maps[i]["atominfofile"] = Path(tempfile.mkstemp(suffix=".dat")[1])

        self._gen_parm7()

        self.inp = Path(tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)[1])
        rmsdfile = Path("rmsd.dat")
        tmp_volumefile = Path(tempfile.mkstemp(suffix=".dat")[1])

        data = {
            "basedir": str(self.basedir),
            "top": self.parm7,
            "traj": self.trajectory,
            "cid": self.probe_id,
            "frame_info": " ".join([str(n) for n in self.frame_info]),
            "ref": self.ref_struct,
            "map_voxel": " ".join([str(n) for n in self.voxel])
            + " gridcenter "
            + " ".join([str(x) for x in box_center]),
            "prefix": self.prefix,
            "maps": maps,
            "rmsdfile": str(rmsdfile),
            "tmp_volumefile": str(tmp_volumefile),
        }

        env = jinja2.Environment(loader=jinja2.FileSystemLoader(f"{os.path.dirname(__file__)}/template"))
        template = env.get_template("cpptraj_pmap.in")
        with open(self.inp, "w") as fout:
            fout.write(template.render(data))
            logger.info(template.render(data))
        command = Command(f"{self.exe} < {self.inp}")
        logger.debug(command)
        try:
            logger.info(command.run())
        except Exception as e:
            Command(f"cat {self.inp}")
            raise e

        for i in range(len(maps)):
            maps[i]["grid"] = self.basedir / f"{self.prefix}_{maps[i]['suffix']}.dx"

        self.frames = len(open(f"{self.basedir}/{rmsdfile}").readlines()) - 1  # -1 for header line
        self.last_volume = float(open(tmp_volumefile).readlines()[-1].split()[1])
        os.system(f"rm {tmp_volumefile}")
        for i in range(len(maps)):
            maps[i]["num_probe_atoms"] = len(open(maps[i]["atominfofile"]).readlines()) - 1  # -1 for header line
            # os.system(f"rm {maps[i]['atominfofile']}")
            # del self.maps[i]["atominfofile"] # it makes errors with multiprocessing
            logger.debug(f"num_probe_atoms {i} {maps[i]['num_probe_atoms']} {maps[i]['atominfofile']}")
        logger.debug(f"{self.trajectory}")

        self.maps = maps

        return self

    def __del__(self):

        if hasattr(self, "inp"):
            logger.debug(f"Cpptraj.inp: {self.inp}")
            if not self.debug:
                os.remove(self.inp)
        if hasattr(self, "parm7"):
            logger.debug(f"Cpptraj.parm7: {self.parm7}")
            if not self.debug:
                os.remove(self.parm7)
