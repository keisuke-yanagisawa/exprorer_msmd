#!/usr/bin/python3

import argparse
import os 

from joblib import Parallel, delayed
import yaml

from script.utilities import util
VERSION = "0.1.0"



def wrapper(index, name, setting_general, setting_input):
    basedirpath = f'{setting_general["workdir"]}/prep{index}'
    inputdirpath = basedirpath+"/input"
    os.system(f"mkdir -p {inputdirpath}")

    tmp_name = "tmp"
    exe_python  = setting_general["executables"]["python"]
    exe_gromacs = setting_general["executables"]["gromacs"]
    # TODO: args.setting_yaml is globally referred
    os.system(f"""
    {exe_python} script/generate_msmd_system.py \
    -setting-yaml {args.setting_yaml} \
    -oprefix {basedirpath}/{tmp_name}_GMX \
    --seed {index}
    """)
    os.system(f"cp {args.setting_yaml} {inputdirpath}")
    
    os.system(f"""
    {exe_python} script/addvirtatom2top.py \
	-i {basedirpath}/{tmp_name}_GMX.top \
	-o {basedirpath}/{tmp_name}_tmp.top \
	-cname {setting_input["probe"]["cid"]} \
	-ovis {basedirpath}/virtual_repulsion.top
    """)

    os.system(f"""
    {exe_python} script/addvirtatom2gro.py \
	-i {basedirpath}/{tmp_name}_GMX.gro \
	-o {basedirpath}/{setting_general["name"]}.gro \
	-cname {setting_input["probe"]["cid"]} \
	-vname "VIS"
    """)

    # gen position restraint files
    os.system(f"""
    {exe_python} script/add_posredefine2top.py \
	-v -res WAT Na+ Cl- CA MG ZN CU {setting_input["probe"]["cid"]} \
	-target protein \
	-gro {basedirpath}/{setting_general["name"]}.gro \
	-i {basedirpath}/{tmp_name}_tmp.top \
	-o {basedirpath}/{setting_general["name"]}.top 
    """)


    os.system(f"""
    {exe_gromacs} trjconv -s {basedirpath}/{setting_general["name"]}.gro \
    -f {basedirpath}/{setting_general["name"]}.gro \
    -o {basedirpath}/{setting_general["name"]}.pdb <<EOF
    0
    EOF
    """)
        
YAML_DIR_PATH_global = None
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="MSMD preparation")
    parser.add_argument("setting_yaml")

    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()
    YAML_DIR_PATH_global = os.path.dirname(
        util.getabsolutepath(args.setting_yaml) 
    )# global variable

    setting = util.parse_yaml(args.setting_yaml)
    indices = util.expand_index(setting["general"]["iter_index"])
    name = setting["general"]["name"]

    setting_general = setting["general"]
    setting_input   = setting["input"]
    Parallel(n_jobs=-1)(
        delayed(wrapper)(idx, name, setting_general, setting_input) 
        for idx in indices
    )
