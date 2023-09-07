import abc
from ..variable import Path
from ..system import SystemInterface, Trajectory


class SimulationInterface(abc.ABC):
    @abc.abstractmethod
    def _create_mdp(self) -> Path:
        pass

    @abc.abstractmethod
    def run(self, initial: SystemInterface) -> Trajectory:
        pass
