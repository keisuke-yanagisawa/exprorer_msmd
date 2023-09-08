from typing import List, Optional
from .interface import SimulationInterface
from ..system import SystemInterface, Trajectory
from ..variable import Path


class SimulationSequence:
    seq: List[SimulationInterface]

    def __init__(self, seq: List[SimulationInterface] = []):
        self.seq = seq

    def run(self, initial: SystemInterface, outdir: Optional[Path] = None) -> Trajectory:
        raise NotImplementedError()


"""
利用例
sys: Final[System] = MSMDSystemCreator.create(CONFIG.INPUT_CONFIG)
sys.add_position_restraint(CONFIG.INPUT_CONFIG.PROBE)
sys.add_pseudo_repulsion()

# うーん。CONFIGをそのままrunできるのは気持ち悪いな
traj: Final[Trajectory] = SIMULATION_CONFIG.SEQUENCE.run(sys)
"""
