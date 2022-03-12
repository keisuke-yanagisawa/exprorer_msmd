import random
import string
import os

import yaml

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
        print("[%s]" % section)
        for option in dat[section]:
            print("  %s : %s" % (option, dat[section][option]))


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

def parse_yaml(yamlpath):
    YAML_PATH = getabsolutepath(yamlpath) 
    YAML_DIR_PATH = os.path.dirname(YAML_PATH)
    with open(YAML_PATH) as fin:
        setting = yaml.safe_load(fin)
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