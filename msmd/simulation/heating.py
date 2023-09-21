from msmd.executable.gromacs import Gromacs
from msmd.jinja2 import render_file
from . import interface
from ..system import System, SystemInterface, Trajectory
from ..variable import Path, Name
from ..unit import Kelvin, Bar, PicoSecond
from .simulation_parameter import NumStep, PressureCoupling
from typing import Final, Dict, Any, Optional


class HeatingStep(interface.SimulationInterface):
    @staticmethod
    def __fill_default(step_config: Dict[str, Any]) -> Dict[str, Any]:
        ret_map = step_config
        ret_map["initial_temp"] = ret_map.get("initial_temp", 0)
        ret_map["target_temp"] = ret_map.get("target_temp", ret_map["temperature"])
        return ret_map

    def _create_mdp(self, seed=-1) -> Path:
        mdp_path = render_file("heating.mdp", suffix=".mdp",
                               define=self.DEFINE,
                               dt=self.DT.get(),
                               nsteps=self.STEPS.get(),
                               nstxtcout=self.NSTXTCOUT.get(),
                               nstlog=self.NSTLOG.get(),
                               nstenergy=self.NSTENERGY.get(),
                               pbc=self.PBC,
                               target_temp=self.TARGET_TEMP.get(),
                               pcoupl=self.PCOUPL.get(),
                               pressure=self.PRESSURE.get(),
                               seed=seed,
                               initial_temp=self.INITIAL_TEMP.get(),
                               duration=self.STEPS.get() * self.DT.get(),
                               )

        return mdp_path

    def __init__(self, step_config: Dict[str, Any]):
        step_config = self.__fill_default(step_config)

        self.NAME: Final[Name] = Name(step_config["name"])
        self.DT: Final[PicoSecond] = PicoSecond(step_config["dt"])
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
        self.PBC: Final[str] = step_config["pbc"]

    def run(self, initial: SystemInterface, outdir: Optional[Path] = None) -> Trajectory:
        input_mdp = self._create_mdp()
        system: System = initial.get_system()
        gromacs = Gromacs()
        traj = gromacs.run(self.NAME, system, input_mdp)
        if outdir is not None:
            gromacs.save(outdir, self.NAME)
            traj.save(outdir, self.NAME)
        return traj

    @property
    def name(self) -> Name:
        return self.NAME
