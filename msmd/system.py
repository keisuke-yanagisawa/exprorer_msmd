import abc
from typing import Final
from .variable import Path, PDBString
from .unit import Angstrom
from .probe.parameter import Probe


class SystemInterface(abc.ABC):
    @abc.abstractmethod
    def get_system(self) -> "System":
        pass

    @abc.abstractmethod
    def save(self, directory: Path) -> None:
        """ファイルに書き出して保存する。"""
        pass


class System(SystemInterface):
    def get_system(self) -> "System":
        return self

    def save(self, directory: Path) -> None:
        pass

    def add_position_restraint(self, probe: Probe) -> None:
        except_residues: Final[list[str]] = [
            "WAT",
            "Na+", "Cl-", "CA", "MG", "ZN", "CU",
            probe.cid
        ]
        raise NotImplementedError()

    def add_pseudo_repulsion(self) -> None:
        sigma: Final[Angstrom] = 20  # 2 nm
        epsilon: Final[float] = 4.184e-6
        # TODO: どのようなルールに基づいて処理をするのか？
        raise NotImplementedError()


class Trajectory(SystemInterface):
    def get_system(self) -> System:
        raise NotImplementedError()

    def save(self, directory: Path) -> None:
        pass


class FrcmodCreator:
    @staticmethod
    def create(probe_config: Probe) -> Path:
        raise NotImplementedError()


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


"""
利用例
sys: Final[System] = MSMDSystemCreator.create(CONFIG.INPUT_CONFIG)
sys.add_position_restraint(CONFIG.INPUT_CONFIG.PROBE)
sys.add_pseudo_repulsion()
"""
