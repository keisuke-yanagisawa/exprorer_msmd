#! /usr/bin/python3


import os
from pathlib import Path
from typing import Dict, List, Optional, Union, cast

import jinja2

from .utilities.logger import logger

VERSION = "1.0.0"


def gen_mdp(protocol_dict: Dict[str, Union[str, float, int]], MD_DIR: Path) -> None:
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
    template = env.get_template(f"./template/{protocol_dict['type']}.mdp")

    if protocol_dict["type"] == "heating":
        if "target_temp" not in protocol_dict:
            protocol_dict["target_temp"] = protocol_dict["temperature"]
        if "initial_temp" not in protocol_dict:
            protocol_dict["initial_temp"] = 0
        nsteps = float(protocol_dict["nsteps"])
        dt = float(protocol_dict["dt"])
        protocol_dict["duration"] = nsteps * dt

    with open(MD_DIR / f"{protocol_dict['name']}.mdp", "w") as fout:
        fout.write(template.render(protocol_dict))


def gen_mdrun_job(
    step_names: List[str],
    name: str,
    path: Path,
    top: Path,
    gro: Path,
    out_traj: Path,
    post_comm: str = ""
) -> None:
    data = {
        "NAME": name,
        "TOP": top,
        "GRO": gro,
        "OUT_TRAJ": out_traj,
        "POST_COMMAND": post_comm,
        "STEP_NAMES": " ".join(step_names),
    }

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
    template = env.get_template("./template/mdrun.sh")
    with open(path, "w") as fout:
        fout.write(template.render(data))
    logger.debug(f"generate {path}")


def prepare_sequence(
    sequence: List[Dict[str, Union[str, float, int]]],
    general: Dict[str, Union[str, float, int]]
) -> List[Dict[str, Union[str, float, int]]]:
    """Prepare MD sequence by combining general and step-specific settings.
    
    Args:
        sequence: List of step-specific settings
        general: General settings to apply to all steps
        
    Returns:
        List of combined settings for each MD step
    """
    ret = []
    for i, step in enumerate(sequence):
        tmp: Dict[str, Union[str, float, int]] = {"define": "", "name": f"step{i+1}"}
        tmp.update(cast(Dict[str, Union[str, float, int]], general))
        tmp.update(cast(Dict[str, Union[str, float, int]], step))
        step = tmp
        ret.append(step)
        logger.info(ret)
    return ret


def prepare_md_files(
    index: int,
    sequence: List[Dict[str, Union[str, float, int]]],
    targetdir: Path,
    jobname: str,
    top: Path,
    gro: Path,
    out_traj: Path
) -> None:
    for step in sequence:
        step["seed"] = index
        gen_mdp(step, targetdir)
    step_names = [str(d["name"]) for d in sequence]
    gen_mdrun_job(step_names, jobname, targetdir / "mdrun.sh", top, gro, out_traj)
