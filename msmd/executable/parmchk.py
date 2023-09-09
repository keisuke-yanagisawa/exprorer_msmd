import tempfile
from typing import Final

from ..probe import parameter
from .command import Executable, Command
from ..variable import Path
from ..standard_library.logging.logger import logger


class Parmchk(object):
    at_indices = {"gaff": 1, "gaff2": 2}

    def __init__(self, exe: Executable, debug=False):
        self.__exe: Final[Executable] = exe
        self.__debug: Final[bool] = debug

    def set(self, probe: "parameter.Probe"):
        self.__probe = probe
        return self

    def run(self) -> Path:
        frcmod = Path(tempfile.mkstemp(suffix=".frcmod")[1])

        command = Command(f"{self.__exe.get()} -i {self.__probe.mol2path.get()} -f mol2 -o {frcmod.path} -s {self.__probe.atomtype.atomtype_index}")
        logger.debug(command)
        logger.info(command.run())
        return frcmod
