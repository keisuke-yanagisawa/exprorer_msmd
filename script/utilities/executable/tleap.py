import os
import tempfile
from pathlib import Path

import jinja2

from .. import const
from ..logger import logger
from .execute import Command


class TLeap(object):
    def __init__(self, debug=False):
        self.exe = os.getenv("TLEAP", "tleap")
        self.debug = debug

    def set(self, cid, probe_path: Path, frcmod: Path, box_path: Path, size, ssbonds, at):
        self.cid = cid
        self.probe_path = probe_path
        self.frcmod = frcmod
        self.box_path = box_path
        self.size = size
        self.ssbonds = ssbonds
        self.at = at
        return self

    def run(self, oprefix: str):
        self.oprefix = oprefix
        self.parm7 = Path(self.oprefix + ".parm7")
        self.rst7 = Path(self.oprefix + ".rst7")

        _, inputfile = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)
        data = {
            "LIGAND_PARAM": f"leaprc.{self.at}",
            "SS_BONDS": self.ssbonds,
            "PROBE_ID": self.cid,
            "PROBE_PATH": str(self.probe_path),
            "OUTPUT": self.oprefix,
            "SYSTEM_PATH": str(self.box_path),
            "PROBE_FRCMOD": str(self.frcmod),
            "SIZE": self.size,
        }

        env = jinja2.Environment(loader=jinja2.FileSystemLoader(f"{os.path.dirname(__file__)}/template"))
        template = env.get_template("leap.in")
        with open(inputfile, "w") as fout:
            fout.write(template.render(data))

        command = Command(f"{self.exe} -f {inputfile}")
        logger.debug(command)
        output = command.run()
        logger.info(output)
        try:
            final_charge_info = [
                s.strip() for s in output.split("\n") if s.strip().startswith("Total unperturbed charge")
            ][0]
        except IndexError as e:
            logger.error(e)
            logger.error(output)
            logger.error(f"cat {self.box_path}")
            logger.error(os.system(f"cat {self.box_path}"))
            return None
        self._final_charge_value = int(float(final_charge_info.split()[-1]))
        return self
