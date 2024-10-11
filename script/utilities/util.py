import collections.abc
import os
import random
import string
from pathlib import Path
from typing import Union

import yaml

from ..setting import (
    Executables,
    GeneralSetting,
    InputSetting,
    MapCreationSetting,
    ProbeProfileSetting,
    ProbeSetting,
    ProfileParameter,
    ProteinSetting,
)
from .logger import logger


def update_dict(d: dict, u: dict):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update_dict(d.get(k, {}), v)  # type: ignore
        else:
            d[k] = v
    return d


def getabsolutepath(path: Path) -> Path:
    """
    Get absolute path from relative path
    """
    path = expandpath(path)
    if not path.is_absolute():  # relative path
        path = os.getcwd() / path
    return path


def expandpath(path: Path) -> Path:
    """
    Expand ~ and $HOME and other environment variables
    """
    path = path.expanduser()
    return Path(os.path.expandvars(path))


def randomname(n: int):
    """
    Generate random string of length n
    """
    randlst = [random.choice(string.ascii_letters + string.digits) for i in range(n)]
    return "".join(randlst)


def dat_dump(dat):
    for section in dat.sections():
        logger.debug("[%s]" % section)
        for option in dat[section]:
            logger.debug("  %s : %s" % (option, dat[section][option]))


# Notation 1-2 => 1,2  5-9:2 => 5,7,9
def expand_index(ind_info: Union[str, int]):
    # single element will be parsed as integer in yaml.safe_load
    # Thus this process is needed
    if isinstance(ind_info, int):
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
            ret.extend(list(range(st, ed + 1)))
        else:
            window, offset = elem.split(":")
            st, ed = window.split("-")
            st, ed, offset = int(st), int(ed), int(offset)
            ret.extend(list(range(st, ed + 1, offset)))
    return ret


def set_default(setting: dict) -> None:
    """
    Set default values (in place) for some fields in setting
    """
    if "multiprocessing" not in setting["general"]:
        setting["general"]["multiprocessing"] = -1
    if "valid_dist" not in setting["map"]:
        setting["map"]["valid_dist"] = 5
    if "threshold" not in setting["probe_profile"]:
        setting["probe_profile"]["threshold"] = 0.001
    if "env_dist" not in setting["probe_profile"]:
        setting["probe_profile"]["resenv"]["env_dist"] = 4

    if "dt" not in setting["exprorer_msmd"]["general"]:
        setting["exprorer_msmd"]["general"]["dt"] = 0.002
    if "temperature" not in setting["exprorer_msmd"]["general"]:
        setting["exprorer_msmd"]["general"]["temperature"] = 300
    if "pressure" not in setting["exprorer_msmd"]["general"]:
        setting["exprorer_msmd"]["general"]["pressure"] = 1.0
    if "num_process_per_gpu" not in setting["general"]:
        setting["general"]["num_process_per_gpu"] = 1


def ensure_compatibility_v1_1(setting: dict):
    """
    Ensure compatibility with exprorer_msmd v1.1
    """

    if "map" not in setting:
        setting["map"] = setting["exprorer_msmd"]["pmap"]
        setting["map"]["snapshot"] = setting["exprorer_msmd"]["pmap"]["snapshots"]
    setting["map"]["snapshot"] = setting["map"]["snapshot"].split("|")[-1]
    if "maps" not in setting["map"]:
        setting["map"]["maps"] = [{"suffix": "nVH", "selector": "(!@VIS)&(!@H*)"}]
    if "map_size" not in setting["map"]:
        setting["map"]["map_size"] = 80
    if "aggregation" not in setting["map"]:
        setting["map"]["aggregation"] = "max"
    if "normalization" not in setting["map"]:
        setting["map"]["normalization"] = "total"


def convert_to_namedtuple(setting: dict) -> None:
    setting["general"]["executables"] = Executables(
        Path(setting["general"]["executables"]["python"]),
        Path(setting["general"]["executables"]["gromacs"]),
        Path(setting["general"]["executables"]["packmol"]),
        Path(setting["general"]["executables"]["tleap"]),
        Path(setting["general"]["executables"]["cpptraj"]),
    )
    setting["general"] = GeneralSetting(
        setting["general"]["iter_index"],
        Path(setting["general"]["workdir"]),
        setting["general"]["name"],
        setting["general"]["executables"],
        setting["general"]["multiprocessing"],
        setting["general"]["num_process_per_gpu"],
    )
    setting["input"]["probe"] = ProbeSetting(
        setting["input"]["probe"]["cid"],
        Path(setting["input"]["probe"]["mol2"]),
        Path(setting["input"]["probe"]["pdb"]),
        setting["input"]["probe"]["atomtype"],
        setting["input"]["probe"]["molar"],
    )
    setting["input"]["protein"] = ProteinSetting(
        Path(setting["input"]["protein"]["pdb"]),
        tuple(setting["input"]["protein"]["ssbond"]),
    )
    setting["input"] = InputSetting(setting["input"]["protein"], setting["input"]["probe"])

    # TODO: ExporerMSMDSetting

    setting["map"] = MapCreationSetting(
        setting["map"]["type"],
        setting["map"]["snapshot"],
        setting["map"]["maps"],
        setting["map"]["map_size"],
        setting["map"]["normalization"],
        setting["map"]["valid_dist"],
        setting["map"]["aggregation"],
    )

    setting["probe_profile"]["profile"]["types"] = [
        ProfileParameter(profile["name"], profile["atoms"]) for profile in setting["probe_profile"]["profile"]["types"]
    ]
    setting["probe_profile"] = ProbeProfileSetting(
        setting["probe_profile"]["resenv"]["map"],
        setting["probe_profile"]["resenv"]["snapshots"],
        setting["probe_profile"]["resenv"]["threshold"],
        setting["probe_profile"]["resenv"]["env_dist"],
        setting["probe_profile"]["resenv"]["align"],
        setting["probe_profile"]["profile"]["types"],
    )


def parse_yaml(yamlpath: Path) -> dict:
    YAML_PATH = getabsolutepath(yamlpath)
    YAML_DIR_PATH = YAML_PATH.parent
    if not YAML_PATH.exists():
        raise FileNotFoundError("YAML file not found: %s" % YAML_PATH)
    if YAML_PATH.is_dir():
        raise IsADirectoryError("Given YAML file path %s is a directory" % YAML_PATH)
    if YAML_PATH.suffix != ".yaml":
        raise ValueError("YAML file must have .yaml extension: %s" % YAML_PATH)
    with YAML_PATH.open() as fin:
        setting: dict = yaml.safe_load(fin)  # type: ignore

    ensure_compatibility_v1_1(setting)
    set_default(setting)

    if "mol2" not in setting["input"]["probe"] or setting["input"]["probe"]["mol2"] is None:
        setting["input"]["probe"]["mol2"] = setting["input"]["probe"]["cid"] + ".mol2"
    if "pdb" not in setting["input"]["probe"] or setting["input"]["probe"]["pdb"] is None:
        setting["input"]["probe"]["pdb"] = setting["input"]["probe"]["cid"] + ".pdb"

    setting["general"]["workdir"] = str(expandpath(Path(setting["general"]["workdir"])))
    setting["general"]["workdir"] = (
        setting["general"]["workdir"]
        if setting["general"]["workdir"].startswith("/")
        or setting["general"]["workdir"].startswith("$HOME")
        or setting["general"]["workdir"].startswith("~")
        else YAML_DIR_PATH / setting["general"]["workdir"]
    )

    setting["input"]["protein"]["pdb"] = str(expandpath(Path(setting["input"]["protein"]["pdb"])))
    setting["input"]["protein"]["pdb"] = (
        setting["input"]["protein"]["pdb"]
        if setting["input"]["protein"]["pdb"].startswith("/")
        or setting["input"]["protein"]["pdb"].startswith("$HOME")
        or setting["input"]["protein"]["pdb"].startswith("~")
        else YAML_DIR_PATH / setting["input"]["protein"]["pdb"]
    )

    setting["input"]["probe"]["pdb"] = str(expandpath(Path(setting["input"]["probe"]["pdb"])))
    setting["input"]["probe"]["pdb"] = (
        setting["input"]["probe"]["pdb"]
        if setting["input"]["probe"]["pdb"].startswith("/")
        or setting["input"]["probe"]["pdb"].startswith("$HOME")
        or setting["input"]["probe"]["pdb"].startswith("~")
        else YAML_DIR_PATH / setting["input"]["probe"]["pdb"]
    )

    setting["input"]["probe"]["mol2"] = str(expandpath(Path(setting["input"]["probe"]["mol2"])))
    setting["input"]["probe"]["mol2"] = (
        setting["input"]["probe"]["mol2"]
        if setting["input"]["probe"]["mol2"].startswith("/")
        or setting["input"]["probe"]["mol2"].startswith("$HOME")
        or setting["input"]["probe"]["mol2"].startswith("~")
        else YAML_DIR_PATH / setting["input"]["probe"]["mol2"]
    )

    if setting["input"]["protein"]["ssbond"] is None:
        setting["input"]["protein"]["ssbond"] = []

    if not "map" in setting["probe_profile"]["resenv"]:
        setting["probe_profile"]["resenv"]["map"] = "nVH"

    convert_to_namedtuple(setting)

    return setting
