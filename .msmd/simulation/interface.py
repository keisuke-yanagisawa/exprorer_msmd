import abc


class SimulationInterface(abc.ABC):
    @abc.abstractmethod
    def run(self) -> None:
        raise NotImplementedError()


class SimulationSettingInterface(abc.ABC):
    pass
