#! /usr/bin/python3


import os
from pathlib import Path
from typing import List

import jinja2

from .utilities.logger import logger

VERSION = "1.0.0"

def gen_mdp(protocol_dict: dict, MD_DIR: Path):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
    if not protocol_dict['type'] in ["minimization", "heating", "equilibration", "production"]:
        raise ValueError(f"Invalid simulation type: {protocol_dict['type']}")
    template = env.get_template(f"./template/{protocol_dict['type']}.mdp")

    if protocol_dict["type"] == "heating":
        if "target_temp" not in protocol_dict:
            protocol_dict["target_temp"] = protocol_dict["temperature"]
        if "initial_temp" not in protocol_dict:
            protocol_dict["initial_temp"] = 0
        protocol_dict["duration"] = protocol_dict["nsteps"] * protocol_dict["dt"]

    with open(MD_DIR / f"{protocol_dict['name']}.mdp", "w") as fout:
        fout.write(template.render(protocol_dict))


def gen_mdrun_job(
    step_names: List[str], name: str, path: Path, top: Path, gro: Path, out_traj: Path, post_comm: str = ""
):
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


def prepare_sequence(sequence, general):
    ret = []
    for i, step in enumerate(sequence):
        tmp = {"define": "", "name": f"step{i+1}"}
        tmp.update(general)
        tmp.update(step)
        step = tmp
        ret.append(step)
        logger.info(ret)
    return ret


def prepare_md_files(
    index: int, sequence: List[dict], targetdir: Path, jobname: str, top: Path, gro: Path, out_traj: Path
):
    for step in sequence:
        step["seed"] = index
        gen_mdp(step, targetdir)
    gen_mdrun_job([d["name"] for d in sequence], jobname, targetdir / "mdrun.sh", top, gro, out_traj)

def run_md_sequence(gpuid: int, simdirpath: Path, exe_gromacs: Path, ncpus: int, jobname: str) -> Path:
    """
    run a simulation sequence with os.system
    """

    # execute simulation
    os.system(f"""
    unset OMP_NUM_THREADS ; \
    export CUDA_VISIBLE_DEVICES="{gpuid}" ; \
    cd {simdirpath} && \
    GMX={exe_gromacs} bash mdrun.sh {ncpus}
    """)

    return simdirpath / f"{jobname}.xtc"