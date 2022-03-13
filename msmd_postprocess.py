#!/usr/bin/python3

import argparse
import os 

import GPUtil
from joblib import Parallel, delayed
import yaml

from script.utilities import util
VERSION = "0.1.0"

from script.utilities.logger import logger

def wrapper(index, name, setting):
    prepdirpath = f'{setting["general"]["workdir"]}/prep{index}'
    sysdirpath = f'{setting["general"]["workdir"]}/system{index}'
    
    if not os.path.exists(prepdirpath):
        raise FileNotFoundError(f"No prep directory: '{prepdirpath}'")
    if not os.path.exists(sysdirpath):
        raise FileNotFoundError(f"No system directory: '{sysdirpath}'")

    tmp_name = "tmp"
    exe_python  = setting["general"]["executables"]["python"]
    exe_gromacs = setting["general"]["executables"]["gromacs"]

    os.system(f"""
    {exe_python} script/genpmap_main.py \
        -basedir {sysdirpath} \
        {setting['general']['yaml']} --debug
    """)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run MSMD simulation")
    parser.add_argument("setting_yaml")

    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    setting = util.parse_yaml(args.setting_yaml)
    indices = util.expand_index(setting["general"]["iter_index"])
    name = setting["general"]["name"]

    Parallel(n_jobs=1)(
        delayed(wrapper)(idx, name, setting) 
        for idx in indices
    )
