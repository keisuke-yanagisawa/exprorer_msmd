from typing import List
from interface import SimulationInterface


class SimulationSequence:
    seq: List[SimulationInterface]

    def __init__(self, seq: List[SimulationInterface] = []):
        self.seq = seq
