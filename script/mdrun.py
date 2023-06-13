#! /usr/bin/python3


from typing import List
import jinja2

import os

from .utilities.logger import logger

VERSION = "1.0.0"


def gen_mdp(protocol_dict: dict,
            MD_DIR: str):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
    template = env.get_template(f"./template/{protocol_dict['type']}.mdp")

    if protocol_dict["type"] == "heating":
        if "target_temp" not in protocol_dict:
            protocol_dict["target_temp"] = protocol_dict["temperature"]
        if "initial_temp" not in protocol_dict:
            protocol_dict["initial_temp"] = 0
        protocol_dict["duration"] = protocol_dict["nsteps"] * protocol_dict["dt"]

    with open(f"{MD_DIR}/{protocol_dict['name']}.mdp", "w") as fout:
        fout.write(template.render(protocol_dict))


def gen_mdrun_job(step_names: List[str],
                  name: str,
                  path: str,
                  top: str,
                  gro: str,
                  out_traj: str,
                  post_comm: str = ""):
    data = {
        "NAME": name,
        "TOP": top,
        "GRO": gro,
        "OUT_TRAJ": out_traj,
        "POST_COMMAND": post_comm,
        "STEP_NAMES": " ".join(step_names)
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


def prepare_md_files(index: int,
                     sequence: List[dict],
                     targetdir: str,
                     jobname: str,
                     top: str,
                     gro: str,
                     out_traj: str):
    for step in sequence:
        step["seed"] = index
        gen_mdp(step, targetdir)
    gen_mdrun_job([d["name"] for d in sequence],
                  jobname, f"{targetdir}/mdrun.sh", top, gro, out_traj)
