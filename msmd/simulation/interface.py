import abc
from ..variable import Path
from ..system import System, Trajectory


class SimulationInterface(abc.ABC):
    @abc.abstractmethod
    def _create_mdp(self) -> Path:
        pass

    @abc.abstractmethod
    def run(self, initial: Trajectory) -> Trajectory:
        pass

    @abc.abstractmethod
    def run_from_system(self, initial: System) -> Trajectory:
        pass
