from typing import List, Optional
from .interface import SimulationInterface
from ..system import System, Trajectory
from ..variable import Path


class SimulationSequence:
    seq: List[SimulationInterface]

    def __init__(self, seq: List[SimulationInterface] = []):
        self.seq = seq

    def run(self, initial: System, outdir: Optional[Path] = None) -> Trajectory:

        first_sim = self.seq[0]
        traj: Trajectory = first_sim.run(initial, outdir)

        other_sims = self.seq[1:]
        for sim in other_sims:
            traj = sim.run(traj, outdir)

        return traj

    def __iter__(self):
        return self.seq.__iter__()
