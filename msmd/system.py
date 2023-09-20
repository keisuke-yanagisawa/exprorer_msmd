import abc
import tempfile
from typing import Final, Optional
from .variable import Path, PDBString
from .unit import Angstrom
from .probe.parameter import Probe
from .parmed import convert as parmed_convert
import parmed as pmd


class SystemInterface(abc.ABC):
    @abc.abstractmethod
    def get_system(self) -> "System":
        pass

    @abc.abstractmethod
    def save(self, directory: Path) -> None:
        """ファイルに書き出して保存する。"""
        pass


class System(SystemInterface):
    @staticmethod
    def create_system_from_strings(top_str: str, gro_str: str) -> "System":
        top: Path = Path(tempfile.mkstemp(suffix=".top")[1])
        gro: Path = Path(tempfile.mkstemp(suffix=".gro")[1])
        with open(top.get(), "w") as f:
            f.write(top_str)
        with open(gro.get(), "w") as f:
            f.write(gro_str)
        return System(top, gro)

    @staticmethod
    def __validation(top: Path, gro: Path) -> None:
        pass
        # 文字列操作で追加する情報があるので、parmedだけでは評価しきれない
        # どのようにvalidationを行うか検討すべき。

    def get_system(self) -> "System":
        return self

    def save(self, dir: Path) -> None:
        pass

    def add_pseudo_repulsion(self, sigma: Angstrom = Angstrom(20), epsilon: float = 4.184e-6) -> None:
        # TODO: どのようなルールに基づいて処理をするのか？
        raise NotImplementedError()

    def __init__(self, top: Path, gro: Path):
        self.__top: Final[Path] = top
        self.__gro: Final[Path] = gro

        self.__validation(self.__top, self.__gro)

    @property
    def top(self) -> Path:
        return self.__top

    @property
    def gro(self) -> Path:
        return self.__gro


class Trajectory(SystemInterface):
    @staticmethod
    def __validation(top: Path, gro: Path, trj: Path) -> None:
        # define as staticmethod to realize const method
        pass

    def get_system(self) -> System:
        raise NotImplementedError()

    def save(self, prefix: Path) -> None:
        pass

    def __init__(self, top: Path, gro: Path, trj: Path):

        self.__top: Final[Path] = top
        self.__gro: Final[Path] = gro
        self.__trj: Final[Path] = trj

    @property
    def top(self) -> Path:
        return self.__top

    @property
    def gro(self) -> Path:
        return self.__gro

    @property
    def trj(self) -> Path:
        return self.__trj


class PDBPreparator:
    @staticmethod
    def read_pdb_file(pdb_path: Path) -> PDBString:
        with open(pdb_path.get(), "r") as f:
            return PDBString(f.read())

    @staticmethod
    def trim_anisou(pdbstring: PDBString) -> PDBString:
        tmp = pdbstring.get().split("\n")
        tmp = [line for line in tmp if line.find("ANISOU") == -1]  # remove ANISOU lines
        return PDBString("\n".join(tmp))

    @staticmethod
    def trim_oxt(pdbstring: PDBString) -> PDBString:
        tmp = pdbstring.get().split("\n")
        tmp = [line for line in tmp if line.find("OXT") == -1]  # remove OXT atoms
        return PDBString("\n".join(tmp))

    @classmethod
    def prepare(cls, pdbstring: PDBString) -> PDBString:
        tmp = cls.trim_anisou(pdbstring)
        tmp = cls.trim_oxt(tmp)
        return tmp


class SystemGenerator():
    def from_amber(self, parm7: Path, rst7: Path) -> System:
        top, gro = parmed_convert(parm7, rst7)
        return System(top, gro)

    def from_gromacs(self, top, gro) -> System:
        raise NotImplementedError()


"""
利用例
sys: Final[System] = MSMDSystemCreator.create(CONFIG.INPUT_CONFIG)
sys.add_position_restraint(CONFIG.INPUT_CONFIG.PROBE)
sys.add_pseudo_repulsion()
"""
