import interface
from ..system import SystemInterface, Trajectory
from ..variable import Path, Name
from .simulation_parameter import NumStep
from typing import Final, Dict, Any


class MinimizationStep(interface.SimulationInterface):
    def _create_mdp(self) -> Path:
        # テンポラリディレクトリパスを作り、そこにファイルを作る
        pass

    def __init__(self, step_config: Dict[str, Any]):
        self.NAME: Final[Name] = Name(step_config["name"])
        self.TITLE: Final[str] = step_config["title"]
        self.STEPS: Final[NumStep] = NumStep(step_config["nsteps"])
        self.DEFINE: Final[str] = step_config["define"]
        self.NSTLOG: Final[NumStep] = NumStep(step_config["nstlog"])

    def run(self, initial: SystemInterface) -> Trajectory:
        # 入力された系に対してMDを実行する
        raise NotImplementedError()


# TODO: heating, Equilibration, Productionも作る
