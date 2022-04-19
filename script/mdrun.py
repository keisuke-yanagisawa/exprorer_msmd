#! /usr/bin/python3

import argparse
import configparser
import pathlib
import subprocess
import shutil
import jinja2
import yaml
import os

VERSION = "1.0.0"


def make_gromacs_directories(parent_dir):
    pathlib.Path("%s/top" % parent_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path("%s/simulation" % parent_dir).mkdir(parents=True, exist_ok=True)
    return parent_dir, "%s/top" % parent_dir, "%s/simulation" % parent_dir


def copy_gromacs_files(inputdir, topdir, name):
    shutil.copy("%s/%s.gro" % (inputdir, name), "%s/%s.gro" % (topdir, name))
    shutil.copy("%s/%s.top" % (inputdir, name), "%s/%s.top" % (topdir, name))
    shutil.copy("%s/index.ndx" % inputdir, "%s/index.ndx" % topdir)

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

def gen_mdrun_job(step_names, name, path, post_comm=""):
    data = {
        "NAME"         : name,
        "POST_COMMAND" : post_comm,
        "STEP_NAMES"   : " ".join(step_names)
    }

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
    template = env.get_template("./template/mdrun.sh")
    with open(path, "w") as fout:
        fout.write(template.render(data))


def update_config(data, var_dict):
    # TODO check key validity
    for key, value in var_dict.items():
        key1, key2 = key.split(":")
        if key1 not in data.keys():
            data[key1] = {}
        data[key1][key2] = value


def validate_config(data):
    # TODO check more things
    if("calc_dir" in data["General"]) and ("output_dir" not in data["General"]):
        print("WARNING: variable General:calc_dir is deprecated. Use General:output_dir instead.")
        data["General"]["output_dir"] = data["General"]["calc_dir"]
    assert("output_dir" in data["General"])

### protocol_dict:define must be rewritten to "" if it is None
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run gromacs jobs automatically")
    parser.add_argument("--version", action="version", version=VERSION)

    # https://gist.github.com/vadimkantorov/37518ff88808af840884355c845049ea
    parser.add_argument("-v", default={},
                        action=type('', (argparse.Action, ),
                                    dict(__call__=lambda a, p, n, v, o:
                                         getattr(n, a.dest).update(dict([v.split('=')])))))

    args = parser.parse_args()

    dat = configparser.ConfigParser()
    update_config(dat, args.v)
    validate_config(dat)
    # print(dat["General"]["name"])

    with open(dat["General"]["protocol_yaml"]) as fin:
        yamldata = yaml.safe_load(fin)

    # 1. make directories
    PARENT_DIR, TOP_DIR, MD_DIR = make_gromacs_directories(dat["General"]["output_dir"])

    for i in range(len(yamldata["exprorer_msmd"]["sequence"])):
        step = yamldata["exprorer_msmd"]["sequence"][i]
        if "name" not in step:
            step["name"] = f"step{i+1}"

        if (not "define" in step) or (step["define"] is None):
            step["define"] = ""
        
        step.update(yamldata["exprorer_msmd"]["general"])
        gen_mdp(step, MD_DIR)
        yamldata["exprorer_msmd"]["sequence"][i] = step # update
        

    gen_mdrun_job([d["name"] for d in yamldata["exprorer_msmd"]["sequence"]],
                  yamldata["general"]["name"],
                  "%s/mdrun.sh" % PARENT_DIR)#,
#                  dat["ProductionRun"]["post_comm"])

    # 1.1 copy raw data
    copy_gromacs_files(dat["General"]["input_dir"], TOP_DIR,
                       yamldata["general"]["name"])

