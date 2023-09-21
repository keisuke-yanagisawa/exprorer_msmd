from . import interface
from ..system import System, Trajectory
from ..variable import Path, Name
from ..unit import Kelvin, Bar
from .simulation_parameter import NumStep, PressureCoupling
from typing import Final, Dict, Any


class HeatingStep(interface.SimulationInterface):
    @staticmethod
    def __fill_default(step_config: Dict[str, Any]) -> Dict[str, Any]:
        ret_map = step_config
        ret_map["initial_temp"] = ret_map.get("initial_temp", 0)
        ret_map["target_temp"] = ret_map.get("target_temp", ret_map["temperature"])
        return ret_map

    def _create_mdp(self) -> Path:
        # テンポラリディレクトリパスを作り、そこにファイルを作る
        raise NotImplementedError()

    def __init__(self, step_config: Dict[str, Any]):
        step_config = self.__fill_default(step_config)

        self.NAME: Final[Name] = Name(step_config["name"])
        self.TITLE: Final[str] = step_config["title"]
        self.STEPS: Final[NumStep] = NumStep(step_config["nsteps"])
        self.DEFINE: Final[str] = step_config["define"]
        self.NSTXTCOUT: Final[NumStep] = NumStep(step_config["nstxtcout"])
        self.NSTLOG: Final[NumStep] = NumStep(step_config["nstlog"])
        self.NSTENERGY: Final[NumStep] = NumStep(step_config["nstenergy"])
        self.PCOUPL: Final[PressureCoupling] = PressureCoupling(step_config["pcoupl"])
        self.INITIAL_TEMP: Final[Kelvin] = Kelvin(step_config["initial_temp"])
        self.TARGET_TEMP: Final[Kelvin] = Kelvin(step_config["target_temp"])
        self.PRESSURE: Final[Bar] = Bar(step_config["pressure"])

    def run(self, initial: Trajectory) -> Trajectory:
        # 入力された系に対してMDを実行する
        raise NotImplementedError()

    def run_from_system(self, initlal: System) -> Trajectory:
        raise RuntimeError("Equilibration step cannot be run from System. Initial minimization steps are required.")
