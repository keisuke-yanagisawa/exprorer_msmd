#! /usr/bin/python3


import os
from pathlib import Path
from typing import Dict, List, Optional, Union, cast

import jinja2

from .utilities.logger import logger

VERSION = "1.0.0"


class MDSimulationManager:
    """Manager class for MD simulation setup and execution."""
    
    def __init__(self, template_dir: Optional[Path] = None) -> None:
        """Initialize MD simulation manager.
        
        Args:
            template_dir: Optional directory containing template files
        """
        self._template_dir = template_dir or Path(os.path.dirname(__file__))
        self._env = jinja2.Environment(loader=jinja2.FileSystemLoader(str(self._template_dir)))
        
    def generate_mdp(self, protocol: Dict[str, Union[str, float, int]], output_dir: Path) -> None:
        """Generate MDP file from protocol settings.
        
        Args:
            protocol: Dictionary containing simulation protocol settings
            output_dir: Directory to write MDP file
        """
        template = self._env.get_template(f"./template/{protocol['type']}.mdp")
        
        if protocol["type"] == "heating":
            if "target_temp" not in protocol:
                protocol["target_temp"] = protocol["temperature"]
            if "initial_temp" not in protocol:
                protocol["initial_temp"] = 0
            nsteps = float(protocol["nsteps"])
            dt = float(protocol["dt"])
            protocol["duration"] = nsteps * dt
            
        with open(output_dir / f"{protocol['name']}.mdp", "w") as fout:
            fout.write(template.render(protocol))
            
    def generate_run_script(
        self,
        step_names: List[str],
        name: str,
        path: Path,
        top: Path,
        gro: Path,
        out_traj: Path,
        post_command: str = ""
    ) -> None:
        """Generate run script for MD simulation.
        
        Args:
            step_names: List of simulation step names
            name: Name of simulation
            path: Path to output script
            top: Path to topology file
            gro: Path to structure file
            out_traj: Path to output trajectory
            post_command: Optional command to run after simulation
        """
        data = {
            "NAME": name,
            "TOP": top,
            "GRO": gro,
            "OUT_TRAJ": out_traj,
            "POST_COMMAND": post_command,
            "STEP_NAMES": " ".join(step_names),
        }
        
        template = self._env.get_template("./template/mdrun.sh")
        with open(path, "w") as fout:
            fout.write(template.render(data))
            logger.debug(f"generate {path}")
            
    def prepare_sequence(
        self,
        sequence: List[Dict[str, Union[str, float, int]]],
        general: Dict[str, Union[str, float, int]]
    ) -> List[Dict[str, Union[str, float, int]]]:
        """Prepare simulation sequence by combining with general settings.
        
        Args:
            sequence: List of step-specific settings
            general: General settings to apply to all steps
            
        Returns:
            List of combined settings for each step
        """
        ret = []
        for i, step in enumerate(sequence):
            tmp: Dict[str, Union[str, float, int]] = {"define": "", "name": f"step{i+1}"}
            tmp.update(general)
            tmp.update(step)
            ret.append(tmp)
            logger.info(ret)
        return ret
            
    def prepare_md_files(
        self,
        index: int,
        sequence: List[Dict[str, Union[str, float, int]]],
        target_dir: Path,
        job_name: str,
        top: Path,
        gro: Path,
        out_traj: Path
    ) -> None:
        """Prepare all files needed for MD simulation.
        
        Args:
            index: Random seed index
            sequence: List of simulation steps
            target_dir: Output directory
            job_name: Name of job
            top: Path to topology file
            gro: Path to structure file
            out_traj: Path to output trajectory
        """
        for step in sequence:
            step["seed"] = index
            self.generate_mdp(step, target_dir)
        self.generate_run_script(
            [str(d["name"]) for d in sequence],
            job_name,
            target_dir / "mdrun.sh",
            top,
            gro,
            out_traj
        )


def gen_mdp(protocol_dict: Dict[str, Union[str, float, int]], MD_DIR: Path) -> None:
    """Generate MDP file (wrapper for backward compatibility)."""
    manager = MDSimulationManager()
    manager.generate_mdp(protocol_dict, MD_DIR)


def gen_mdrun_job(
    step_names: List[str],
    name: str,
    path: Path,
    top: Path,
    gro: Path,
    out_traj: Path,
    post_comm: str = ""
) -> None:
    """Generate run script (wrapper for backward compatibility)."""
    manager = MDSimulationManager()
    manager.generate_run_script(step_names, name, path, top, gro, out_traj, post_comm)


def prepare_sequence(
    sequence: List[Dict[str, Union[str, float, int]]],
    general: Dict[str, Union[str, float, int]]
) -> List[Dict[str, Union[str, float, int]]]:
    """Prepare MD sequence (wrapper for backward compatibility)."""
    manager = MDSimulationManager()
    return manager.prepare_sequence(sequence, general)


def prepare_md_files(
    index: int,
    sequence: List[Dict[str, Union[str, float, int]]],
    targetdir: Path,
    jobname: str,
    top: Path,
    gro: Path,
    out_traj: Path
) -> None:
    """Prepare MD files (wrapper for backward compatibility)."""
    manager = MDSimulationManager()
    manager.prepare_md_files(index, sequence, targetdir, jobname, top, gro, out_traj)
