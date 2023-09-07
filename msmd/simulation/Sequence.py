from typing import List
from interface import SimulationInterface
from ..system import SystemInterface, Trajectory


class SimulationSequence:
    seq: List[SimulationInterface]

    def __init__(self, seq: List[SimulationInterface] = []):
        self.seq = seq

    def run(self, initial: SystemInterface) -> Trajectory:
        pass
