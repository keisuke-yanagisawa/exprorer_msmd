from . import interface
from ..system import SystemInterface, Trajectory
from ..variable import Path, Name
from ..unit import Kelvin, Bar
from .simulation_parameter import NumStep, PressureCoupling
from typing import Final, Dict, Any


class EquilibrationStep(interface.SimulationInterface):
    def _create_mdp(self) -> Path:
        # テンポラリディレクトリパスを作り、そこにファイルを作る
        pass

    def __init__(self, step_config: Dict[str, Any]):
        self.NAME: Final[Name] = Name(step_config["name"])
        self.TITLE: Final[str] = step_config["title"]
        self.STEPS: Final[NumStep] = NumStep(step_config["nsteps"])
        self.DEFINE: Final[str] = step_config["define"]
        self.NSTXTCOUT: Final[NumStep] = NumStep(step_config["nstxtcout"])
        self.NSTLOG: Final[NumStep] = NumStep(step_config["nstlog"])
        self.NSTENERGY: Final[NumStep] = NumStep(step_config["nstenergy"])
        self.PCOUPL: Final[PressureCoupling] = PressureCoupling(step_config["pcoupl"])
        self.TEMPERATURE: Final[Kelvin] = Kelvin(step_config["temperature"])
        self.PRESSURE: Final[Bar] = Bar(step_config["pressure"])

    def run(self, initial: SystemInterface) -> Trajectory:
        # 入力された系に対してMDを実行する
        raise NotImplementedError()
