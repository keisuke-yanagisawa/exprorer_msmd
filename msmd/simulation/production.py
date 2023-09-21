from msmd.executable.gromacs import Gromacs
from msmd.jinja2 import render_file
from . import interface
from ..system import System, SystemInterface, Trajectory
from ..variable import Path, Name
from ..unit import Kelvin, Bar, PicoSecond
from .simulation_parameter import NumStep, PressureCoupling
from typing import Final, Dict, Any


class ProductionStep(interface.SimulationInterface):
    def _create_mdp(self, seed=-1) -> Path:
        mdp_path = render_file("production.mdp", suffix=".mdp",
                               define=self.DEFINE,
                               dt=self.DT.get(),
                               nsteps=self.STEPS.get(),
                               nstxtcout=self.NSTXTCOUT.get(),
                               nstenergy=self.NSTENERGY.get(),
                               nstlog=self.NSTLOG.get(),
                               seed=seed,
                               temperature=self.TEMPERATURE.get(),
                               pcoupl=self.PCOUPL.get(),
                               pressure=self.PRESSURE.get(),
                               pbc=self.PBC)

        return mdp_path

    def __init__(self, step_config: Dict[str, Any]):
        self.NAME: Final[Name] = Name(step_config["name"])
        self.DT: Final[PicoSecond] = PicoSecond(step_config["dt"])
        self.STEPS: Final[NumStep] = NumStep(step_config["nsteps"])
        self.DEFINE: Final[str] = step_config["define"]
        self.NSTXTCOUT: Final[NumStep] = NumStep(step_config["nstxtcout"])
        self.NSTLOG: Final[NumStep] = NumStep(step_config["nstlog"])
        self.NSTENERGY: Final[NumStep] = NumStep(step_config["nstenergy"])
        self.PCOUPL: Final[PressureCoupling] = PressureCoupling(step_config["pcoupl"])
        self.TEMPERATURE: Final[Kelvin] = Kelvin(step_config["temperature"])
        self.PRESSURE: Final[Bar] = Bar(step_config["pressure"])
        self.PBC: Final[str] = step_config["pbc"]

        self.TITLE: Final[str] = step_config["title"]

    def run(self, initial: SystemInterface) -> Trajectory:
        input_mdp = self._create_mdp()
        system: System = initial.get_system()
        return Gromacs().run(self.NAME, system, input_mdp)
