#! /usr/bin/python3

import argparse
import jinja2
import yaml
import os

from .utilities.logger import logger

VERSION = "1.0.0"


def gen_mdp(protocol_dict, MD_DIR):
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


def gen_mdrun_job(step_names, name, path, top, gro, out_traj, post_comm=""):
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


def prepare_md_files(sequence, targetdir, jobname, top, gro, out_traj):
    for step in sequence:
        gen_mdp(step, targetdir)
    gen_mdrun_job([d["name"] for d in sequence],
                  jobname, f"{targetdir}/mdrun.sh", top, gro, out_traj)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run gromacs jobs automatically")
    parser.add_argument("--version", action="version", version=VERSION)

    parser.add_argument("-dir")
    parser.add_argument("yaml")
    args = parser.parse_args()

    with open(args.yaml) as fin:
        yamldata = yaml.safe_load(fin)

    yamldata["exprorer_msmd"]["sequence"] = prepare_sequence(
        yamldata["exprorer_msmd"]["sequence"],
        yamldata["exprorer_msmd"]["general"]
    )

    prepare_md_files(
        yamldata["exprorer_msmd"]["sequence"],
        args.dir,
        yamldata["general"]["name"]
    )
