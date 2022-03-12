#!/usr/bin/python3

import argparse
import os 

import GPUtil
from joblib import Parallel, delayed
import yaml

from script.utilities import util
VERSION = "0.1.0"

from script.utilities.logger import logger

def wrapper(index, name, setting, gpuid, ncpus):
    prepdirpath = f'{setting["general"]["workdir"]}/prep{index}'
    sysdirpath = f'{setting["general"]["workdir"]}/system{index}'
    topdirpath = f'{sysdirpath}/top'
    
    if not os.path.exists(prepdirpath):
        raise FileNotFoundError(f"No prep directory: '{prepdirpath}'")
    if not os.path.exists(sysdirpath):
        logger.info(f"make {sysdirpath}")
        os.system(f"mkdir -p {sysdirpath}")
    if not os.path.exists(topdirpath):
        logger.info(f"make {topdirpath}")
        os.system(f"mkdir -p {topdirpath}")

    tmp_name = "tmp"
    exe_python  = setting["general"]["executables"]["python"]
    exe_gromacs = setting["general"]["executables"]["gromacs"]

    # TODO: args.setting_yaml is globally referred
    os.system(f"""
    {exe_python} script/mdrun.py \
	DUMMY.conf \
	script/mdrun.sh \
	-v General:input_dir={prepdirpath} \
	-v General:output_dir={sysdirpath} \
	-v General:name={setting["general"]["name"]} \
    -v General:protocol_yaml={args.setting_yaml}
    # """)

    os.system(f"""
    cd {topdirpath} &&
    {exe_gromacs} make_ndx -f {prepdirpath}/{setting["general"]["name"]}.gro << EOF
    q
    EOF
    """)
    
    os.system(f"""
    unset OMP_NUM_THREADS
    export CUDA_VISIBLE_DEVICES="{gpuid}"
    cd {sysdirpath} && \
    GMX={exe_gromacs} bash mdrun.sh {ncpus}
    """)
        
YAML_DIR_PATH_global = None
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run MSMD simulation")
    parser.add_argument("setting_yaml")

    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()
    YAML_DIR_PATH_global = os.path.dirname(
        util.getabsolutepath(args.setting_yaml) 
    )# global variable

    setting = util.parse_yaml(args.setting_yaml)
    indices = util.expand_index(setting["general"]["index"])
    name = setting["general"]["name"]

    gpuids = GPUtil.getAvailable()
    ngpus = len(gpuids)
    if ngpus == 0: # there is only CPU
        gpuids = [-1]
        ngpus = 1
    ncpus = os.cpu_count()
    ncpus_per_gpu = ncpus // ngpus

    gpuids = (gpuids * len(indices))[:len(indices)]

    Parallel(n_jobs=ngpus)(
        delayed(wrapper)(idx, name, setting, gpuid, ncpus_per_gpu) 
        for idx, gpuid in zip(indices, gpuids)
    )
