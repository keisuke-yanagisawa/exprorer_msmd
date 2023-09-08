from typing import Final, List
from ..variable import Path, Name
from ..config import Protein, _PMAPConfig, _PMAPTypeConfig
from ..system import Trajectory
from gridData import Grid


class PMAP:
    """
    PMAPの情報を保持するクラス
    """

    def __init__(self, grid: Grid, name: Name):
        self.__grid: Final[Grid] = grid
        self.__name: Final[Name] = name
        pass


class PMAPCreator:
    def __create_single_pmap(self, pmap_type: _PMAPTypeConfig) -> None:
        raise NotImplementedError()

    @staticmethod
    def create(trajectory: Trajectory, pmap_config: _PMAPConfig) -> List[PMAP]:
        raise NotImplementedError()

    def __init__(self, outdir: Path, original_protein: Protein):
        self.__outdir: Final[Path] = outdir
        self.__reference_protein: Final[Protein] = original_protein


"""
利用例
sys: Final[System] = MSMDSystemCreator.create(CONFIG.INPUT_CONFIG)
sys.add_position_restraint(CONFIG.INPUT_CONFIG.PROBE)
sys.add_pseudo_repulsion()

# うーん。CONFIGをそのままrunできるのは気持ち悪いな
traj: Final[Trajectory] = SIMULATION_CONFIG.SEQUENCE.run(sys)

pc = PMAPCreator(Path("directory"), CONFIG.INPUT_CONFIG.PROTEIN)
pmap_list = pc.create(traj, CONFIG.PMAP_CONFIG)
"""
