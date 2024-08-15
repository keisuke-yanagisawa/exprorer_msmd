from pathlib import Path
from typing import NamedTuple


class Executables(NamedTuple):
    python: Path = Path("python")
    gromacs: Path = Path("gmx")
    packmol: Path = Path("packmol")
    tleap: Path = Path("tleap")
    cpptraj: Path = Path("cpptraj")


class GeneralSetting(NamedTuple):
    iter_index: list[int]
    workdir: Path
    name: str
    executables: Executables
    multiprocessing: int = -1  # TODO: duplicated?
    num_process_per_gpu: int = 1


class ProteinSetting(NamedTuple):
    pdb: Path
    ssbond: tuple[tuple[int, int], ...] = tuple()


class ProbeSetting(NamedTuple):
    cid: str
    mol2: Path
    pdb: Path
    atomtype: str = "gaff2"
    molar: float = 0.25


class InputSetting(NamedTuple):
    protein: ProteinSetting
    probe: ProbeSetting

