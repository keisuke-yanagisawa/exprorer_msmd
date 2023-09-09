import shutil
import tempfile
from typing import Final

from ..jinja2 import render_str
from ..probe.parameter import Probe
from ..protein.parameter import Protein
from ..unit import Angstrom
from .command import Executable, Command
from ..variable import PDBString, Path
from ..standard_library.logging.logger import logger
import warnings


class Packmol(object):
    @staticmethod
    def __shorten_path(path: Path, ext: str) -> Path:
        ret_path = Path(tempfile.mkstemp(suffix=f".{ext}")[1])
        shutil.copy2(path.get(), ret_path.get())
        return ret_path

    def _validation(self) -> None:
        if self.__num_probes < 0:
            raise ValueError(f"The number of probes: {self.__num_probes} is not positive")
        if self.__num_probes == 0:
            warnings.warn("There is no probe molecule.", RuntimeWarning)

    def __init__(self, exe: Executable, debug=False):
        self.exe: Final[Executable] = exe
        self.debug: Final[bool] = debug

    def set(self, protein: Protein, probe: Probe, system_length: Angstrom, num_probes: int) -> "Packmol":
        self.__protein = protein
        self.__probe = probe
        self.__system_length = system_length
        self.__num_probes = num_probes

        self._validation()
        return self

    def run(self, seed=-1) -> PDBString:
        box_pdb: Final[Path] = Path(tempfile.mkstemp(suffix=".pdb")[1])
        tmp_protein_pdb: Final[Path] = self.__shorten_path(self.__protein.pdbpath, "pdb")
        tmp_probe_pdb: Final[Path] = self.__shorten_path(self.__probe.pdbpath, "pdb")
        self.__seed = seed

        probes = [{
            "pdb": tmp_probe_pdb.path,
            "num": self.__num_probes,
            "size": self.__system_length.get() / 2
        }]

        data = {
            "output": box_pdb.path,
            "prot": tmp_protein_pdb.path,
            "seed": self.__seed,
            "probes": [probe for probe in probes if probe["num"] > 0],
        }

        input_str = render_str("packmol.in", **data)

        _, inp = tempfile.mkstemp(suffix=".in")
        with open(inp, "w") as fout:
            fout.write(input_str)
        if self.debug:
            logfile = "packmol.log"
        else:
            _, logfile = tempfile.mkstemp(suffix=".log")

        command = Command(f"{self.exe.get()} < {inp} > {logfile}")
        logger.debug(command)
        command.run()

        stdout_str = open(logfile).read()
        if stdout_str.find("ENDED WITHOUT PERFECT PACKING") != -1:
            raise RuntimeError(f"packmol failed: \n{stdout_str}")
        logger.info(stdout_str)
        return PDBString(open(box_pdb.path).read())
