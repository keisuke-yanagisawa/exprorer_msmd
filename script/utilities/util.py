import math
import random
import string
import os

import GPUtil
import yaml

from .logger import logger


def getabsolutepath(path):
    path = expandpath(path)
    if not path.startswith("/"): # relative path
        path = os.getcwd() + "/" + path
    return path

def expandpath(path):
    path = os.path.expanduser(path)
    return os.path.expandvars(path)


def randomname(n):
    randlst = [random.choice(string.ascii_letters + string.digits) for i in range(n)]
    return "".join(randlst)


def dat_dump(dat):
    for section in dat.sections():
        logger.debug("[%s]" % section)
        for option in dat[section]:
            logger.debug("  %s : %s" % (option, dat[section][option]))


# Notation 1-2 => 1,2  5-9:2 => 5,7,9
def expand_index(ind_info):
    # single element will be parsed as integer in yaml.safe_load
    # Thus this process is needed
    if type(ind_info) == int:
        return [ind_info]
    
    ret = []
    elems = ind_info.split(",")
    for elem in elems:
        if elem.find("-") == -1:
            elem = int(elem)
            ret.append(elem)
        elif elem.find(":") == -1:
            st, ed = elem.split("-")
            st, ed = int(st), int(ed)
            ret.extend(list(range(st, ed+1)))
        else:
            window, offset = elem.split(":")
            st, ed = window.split("-")
            st, ed, offset = int(st), int(ed), int(offset)
            ret.extend(list(range(st, ed+1, offset)))
    return ret

def get_gpuids():
    gpuids = set(GPUtil.getAvailable(limit=math.inf))
    logger.info(f"{len(gpuids)} GPUs are detected")
    if os.getenv("CUDA_VISIBLE_DEVICES") is not None:
        logger.info(f"CUDA_VISIBLE_DEVICES detected")
        cvd = [int(s) for s in os.getenv("CUDA_VISIBLE_DEVICES", default="").split(",")]
        gpuids &= set(cvd)
        gpuids = list(gpuids)

    if len(gpuids) == 0:
        logger.warn(f"No GPU is allowed/existed to use")
        logger.warn(f"Switch to CPU-only mode, it greatly decreases the simulation speed")
    else:
        logger.info(f"GPU IDs of {gpuids} will be used")
    
    ngpus = len(gpuids)
    if ngpus == 0: # there is only CPU
        gpuids = [-1]
        ngpus = 1
    
    return list(gpuids)

def set_default(setting):
    if not "multiprocessing" in setting["general"]:
        setting["general"]["multiprocessing"] = True    
    if not "valid_dist" in setting["map"]:
        setting["map"]["valid_dist"] = 5
    if not "env_dist" in setting["probe_profile"]:
        setting["probe_profile"]["env_dist"] = 4
    if not "dt" in setting["exprorer_msmd"]["general"]:
        setting["exprorer_msmd"]["general"]["dt"] = 0.002
    if not "temperature" in setting["exprorer_msmd"]["general"]:
        setting["exprorer_msmd"]["general"]["temperature"] = 300
    if not "pressure" in setting["exprorer_msmd"]["general"]:
        setting["exprorer_msmd"]["general"]["pressure"] = 1.0
    if not "seed" in setting["exprorer_msmd"]["general"]:
        setting["exprorer_msmd"]["general"]["seed"] = -1

def ensure_compatibility_v1_1(setting):
    if not "map" in setting:
        setting["map"] = setting["exprorer_msmd"]["pmap"]
        setting["map"]["snapshot"] = setting["exprorer_msmd"]["pmap"]["snapshots"]
    setting["map"]["snapshot"] = setting["map"]["snapshot"].split("|")[-1]
    if not "maps" in setting["map"]:
        setting["map"]["maps"] = [
            {"suffix": "nVH", "selector": "(!@VIS)&(!@H*)"}
        ]
    if not "map_size" in setting["map"]:
        setting["map"]["map_size"] = 80
    if not "aggregation" in setting["map"]:
        setting["map"]["aggregation"] = "max"
    if not "normalization" in setting["map"]:
        setting["map"]["normalization"] = "total"

def parse_yaml(yamlpath):
    YAML_PATH = getabsolutepath(yamlpath) 
    YAML_DIR_PATH = os.path.dirname(YAML_PATH)
    with open(YAML_PATH) as fin:
        setting = yaml.safe_load(fin)

    ensure_compatibility_v1_1(setting)
    set_default(setting)

    if not "mol2" in setting["input"]["probe"] \
       or setting["input"]["probe"]["mol2"] is None:
        setting["input"]["probe"]["mol2"] \
            = setting["input"]["probe"]["cid"]+".mol2"
    if not "pdb" in setting["input"]["probe"] \
       or setting["input"]["probe"]["pdb"] is None:
        setting["input"]["probe"]["pdb"] \
            = setting["input"]["probe"]["cid"]+".pdb"

    setting["general"]["workdir"] = setting["general"]["workdir"] \
        if setting["general"]["workdir"].startswith("/") \
        else YAML_DIR_PATH + "/" + setting["general"]["workdir"]
    setting["input"]["protein"]["pdb"] = setting["input"]["protein"]["pdb"] \
        if setting["input"]["protein"]["pdb"].startswith("/") \
        else YAML_DIR_PATH + "/" + setting["input"]["protein"]["pdb"]
    setting["input"]["probe"]["pdb"] = setting["input"]["probe"]["pdb"] \
        if setting["input"]["probe"]["pdb"].startswith("/") \
        else YAML_DIR_PATH + "/" + setting["input"]["probe"]["pdb"]
    setting["input"]["probe"]["mol2"] = setting["input"]["probe"]["mol2"] \
        if setting["input"]["probe"]["mol2"].startswith("/") \
        else YAML_DIR_PATH + "/" + setting["input"]["probe"]["mol2"]

    if setting["input"]["protein"]["ssbond"] is None:
        setting["input"]["protein"]["ssbond"] = []

    setting["general"]["yaml"] = YAML_PATH

    return setting