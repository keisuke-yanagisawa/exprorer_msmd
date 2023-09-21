from typing import List, Optional
from .interface import SimulationInterface
from ..system import System, Trajectory
from ..variable import Path


class SimulationSequence:
    seq: List[SimulationInterface]

    def __init__(self, seq: List[SimulationInterface] = []):
        self.seq = seq

    def run(self, initial: System, outdir: Optional[Path] = None) -> Trajectory:
        do_output_files = outdir is not None

        first_sim = self.seq[0]
        other_sims = self.seq[1:]

        traj: Trajectory = first_sim.run_from_system(initial)
        if do_output_files:
            traj.save(outdir)

        for sim in other_sims:
            traj = sim.run(traj)
            if do_output_files:
                traj.save(outdir)

        return traj
