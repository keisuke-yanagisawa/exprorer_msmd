import tempfile
import os
from .execute import Command
from .. import const
from ..logger import logger


class Parmchk(object):
    def __init__(self, debug=False):
        self.at_indices = {"gaff": 1, "gaff2": 2}
        self.exe = os.getenv("PARMCHK", "parmchk2")
        self.debug = debug

    def set(self, mol2, at):
        if at not in self.at_indices:
            raise ValueError(f"atomtype {at} is not supported")
        if not os.path.exists(mol2):
            raise FileNotFoundError(f"{mol2} does not exist")
        if not os.path.splitext(mol2)[1][1:] == "mol2":
            raise ValueError(f"{mol2} is not mol2 file")
        self.at_id = self.at_indices[at]
        self.mol2 = mol2
        return self

    def run(self, frcmod=None):
        self.frcmod = frcmod
        if self.frcmod is None:
            _, self.frcmod = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_FRCMOD)

        os.makedirs(os.path.dirname(self.frcmod), exist_ok=True)  # mkdir -p
        command = Command(f"{self.exe} -i {self.mol2} -f mol2 -o {self.frcmod} -s {self.at_id}")
        logger.debug(command)
        logger.info(command.run())
        return self
