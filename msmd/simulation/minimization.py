import tempfile
from . import interface
from ..system import System, SystemInterface, Trajectory
from ..variable import Path, Name
from .simulation_parameter import NumStep
from typing import Final, Dict, Any
from msmd.executable.gromacs import Gromacs
from msmd.jinja2 import render_file


class MinimizationStep(interface.SimulationInterface):
    def _create_mdp(self) -> Path:
        mdp_path = render_file("minimization.mdp", suffix=".mdp",
                               define=self.DEFINE,
                               nsteps=self.STEPS.get(),
                               nstlog=self.NSTLOG.get(),
                               pbc=self.PBC)

        return mdp_path

    def __init__(self, step_config: Dict[str, Any]):
        self.NAME: Final[Name] = Name(step_config["name"])
        self.TITLE: Final[str] = step_config["title"]
        self.STEPS: Final[NumStep] = NumStep(step_config["nsteps"])
        self.DEFINE: Final[str] = step_config["define"]
        self.NSTLOG: Final[NumStep] = NumStep(step_config["nstlog"])
        self.PBC: Final[str] = step_config["pbc"]

    def run(self, initial: SystemInterface) -> Trajectory:
        input_mdp = self._create_mdp()
        system: System = initial.get_system()
        return Gromacs().run(self.NAME, system, input_mdp)
