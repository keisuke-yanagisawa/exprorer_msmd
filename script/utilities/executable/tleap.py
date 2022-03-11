import tempfile
import os
from subprocess import getoutput as gop
import jinja2
from utilities import const
from utilities.executable.execute import Command


class TLeap(object):
    def __init__(self, debug=False):
        self.exe = os.getenv("TLEAP", "tleap")
        self.debug = debug

    def set(self, cids, cosolv_paths, frcmods, box_path, size, ssbonds, at):
        self.cids = cids
        self.cosolv_paths = cosolv_paths
        self.frcmods = frcmods
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
            "COSOLVENTS": [{"ID": cid, "PATH": cpath}
                        for cid, cpath in zip(self.cids, self.cosolv_paths)],
            "OUTPUT": self.oprefix,
            "SYSTEM_PATH": self.box_path,
            "COSOLV_FRCMODS": self.frcmods,
            "SIZE": self.size
        }

        env = jinja2.Environment(loader=jinja2.FileSystemLoader(f"{os.path.dirname(__file__)}/template"))
        template = env.get_template("leap.in")
        with open(inputfile, "w") as fout:
            fout.write(template.render(data))

        command = Command(f"{self.exe} -f {inputfile}")
        if self.debug:
          print(command)
        output = command.run()
        print(output)
        final_charge_info = [s.strip() for s in output.split("\n")
                            if s.strip().startswith("Total unperturbed charge")][0]
        self._final_charge_value = float(final_charge_info.split()[-1])
        return self