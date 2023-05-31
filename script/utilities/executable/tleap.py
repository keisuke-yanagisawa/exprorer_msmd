import tempfile
import os
import jinja2
from .. import const
from .execute import Command
from ..logger import logger


class TLeap(object):
    def __init__(self, debug=False):
        self.exe = os.getenv("TLEAP", "tleap")
        self.debug = debug

    def set(self, cid, probe_path, frcmod, box_path, size, ssbonds, at):
        self.cid = cid
        self.probe_path = probe_path
        self.frcmod = frcmod
        self.box_path = box_path
        self.size = size
        self.ssbonds = ssbonds
        self.at = at
        return self

    def run(self, oprefix):
        self.oprefix = oprefix
        self.parm7 = self.oprefix + ".parm7"
        self.rst7 = self.oprefix + ".rst7"

        _, inputfile = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)
        data = {
            "LIGAND_PARAM": f"leaprc.{self.at}",
            "SS_BONDS": self.ssbonds,
            "PROBE_ID": self.cid,
            "PROBE_PATH": self.probe_path,
            "OUTPUT": self.oprefix,
            "SYSTEM_PATH": self.box_path,
            "PROBE_FRCMOD": self.frcmod,
            "SIZE": self.size
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
            final_charge_info = [s.strip() for s in output.split("\n")
                                 if s.strip().startswith("Total unperturbed charge")][0]
        except IndexError as e:
            logger.error(e)
            logger.error(output)
            logger.error(f"cat {self.box_path}")
            logger.error(os.system(f"cat {self.box_path}"))
            return None
        self._final_charge_value = float(final_charge_info.split()[-1])
        return self
