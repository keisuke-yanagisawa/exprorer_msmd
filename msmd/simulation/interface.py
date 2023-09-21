import abc
from typing import Optional
from ..variable import Name, Path
from ..system import SystemInterface, Trajectory


class SimulationInterface(abc.ABC):
    @abc.abstractmethod
    def _create_mdp(self, seed=-1) -> Path:
        pass

    @abc.abstractmethod
    def run(self, initial: SystemInterface, outdir: Optional[Path] = None) -> Trajectory:
        pass

    @property
    @abc.abstractmethod
    def name(self) -> Name:
        pass
