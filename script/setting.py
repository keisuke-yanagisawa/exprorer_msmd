from pathlib import Path
from typing import Literal, NamedTuple


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


class MapCreationSetting(NamedTuple):
    type: Literal["pmap"]  # TODO: pmap以外もあるはず
    snapshot: str
    maps: list[dict] # TODO: mapの型を準備すべき。ただし、後から情報が付け足されるコードが存在するため要検討
    map_size: int = 80  # angstrom
    normalization: Literal["total", "snapshot", "GFE"] = "total"
    valid_dist: float = 5  # angstrom
    aggregation: Literal["max"] = "max"


class ProfileParameter(NamedTuple):
    name: str
    atoms: tuple[tuple[str, str], ...]


class ProbeProfileSetting(NamedTuple):
    map: str
    snapshots: str
    threshold: float = 0.001  # 0.1 percent
    env_dist: float = 4.0  # angstrom
    align: list[str] = [" C1 ", " C2 ", " C3 ", " O1 "]
    profiles: list[ProfileParameter] = []

