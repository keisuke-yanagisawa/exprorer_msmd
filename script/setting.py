import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import yaml

from .utilities.util import expandpath, getabsolutepath, update_dict


def ensure_compatibility_v1_1(setting: Dict[str, Any]) -> None:
    """
    Ensure compatibility with exprorer_msmd v1.1
    """

    if "pmap" in setting["exprorer_msmd"]:
        update_dict(setting["map"], setting["exprorer_msmd"]["pmap"])  # because of renaming
    setting["map"]["snapshot"] = setting["map"]["snapshot"].split("|")[-1]  # use only one snapshot data


def parse_yaml(yamlpath: Path) -> Dict[str, Any]:
    setting: dict = {
        "general": {
            "workdir": Path(""),
            "multiprocessing": -1,
            "num_process_per_gpu": 1,
        },
        "input": {
            "protein": {
                "pdb": Path(""),
            },
            "probe": {
                "cid": "",
            },
        },
        "exprorer_msmd": {
            "general": {
                "dt": 0.002,
                "temperature": 300,
                "pressure": 1.0,
            },
        },
        "map": {
            "snapshot": "",
            "valid_dist": 5.0,
            "map_size": 80,
            "normalization": "total",
            "aggregation": "max",
            "maps": [
                {
                    "suffix": "nVH",
                    "selector": "(!@VIS)&(!@H*)",
                }
            ],
        },
        "probe_profile": {
            "threshold": 0.001,
            "resenv": {
                "env_dist": 4.0,
            },
        },
    }
    YAML_PATH = getabsolutepath(yamlpath)
    YAML_DIR_PATH = YAML_PATH.parent
    if not YAML_PATH.exists():
        raise FileNotFoundError("YAML file not found: %s" % YAML_PATH)
    if YAML_PATH.is_dir():
        raise IsADirectoryError("Given YAML file path %s is a directory" % YAML_PATH)
    if not os.path.splitext(YAML_PATH)[1][1:] == "yaml":
        raise ValueError("YAML file must have .yaml extension: %s" % YAML_PATH)
    with YAML_PATH.open() as fin:
        yaml_dict: dict = yaml.safe_load(fin)  # type: ignore
        if yaml_dict is not None:
            update_dict(setting, yaml_dict)

    ensure_compatibility_v1_1(setting)

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
    setting["general"]["workdir"] = Path(setting["general"]["workdir"])

    setting["input"]["protein"]["pdb"] = str(expandpath(Path(setting["input"]["protein"]["pdb"])))
    setting["input"]["protein"]["pdb"] = (
        setting["input"]["protein"]["pdb"]
        if setting["input"]["protein"]["pdb"].startswith("/")
        or setting["input"]["protein"]["pdb"].startswith("$HOME")
        or setting["input"]["protein"]["pdb"].startswith("~")
        else YAML_DIR_PATH / setting["input"]["protein"]["pdb"]
    )
    setting["input"]["protein"]["pdb"] = Path(setting["input"]["protein"]["pdb"])

    setting["input"]["probe"]["pdb"] = str(expandpath(Path(setting["input"]["probe"]["pdb"])))
    setting["input"]["probe"]["pdb"] = (
        setting["input"]["probe"]["pdb"]
        if setting["input"]["probe"]["pdb"].startswith("/")
        or setting["input"]["probe"]["pdb"].startswith("$HOME")
        or setting["input"]["probe"]["pdb"].startswith("~")
        else YAML_DIR_PATH / setting["input"]["probe"]["pdb"]
    )
    setting["input"]["probe"]["pdb"] = Path(setting["input"]["probe"]["pdb"])

    setting["input"]["probe"]["mol2"] = str(expandpath(Path(setting["input"]["probe"]["mol2"])))
    setting["input"]["probe"]["mol2"] = (
        setting["input"]["probe"]["mol2"]
        if setting["input"]["probe"]["mol2"].startswith("/")
        or setting["input"]["probe"]["mol2"].startswith("$HOME")
        or setting["input"]["probe"]["mol2"].startswith("~")
        else YAML_DIR_PATH / setting["input"]["probe"]["mol2"]
    )
    setting["input"]["probe"]["mol2"] = Path(setting["input"]["probe"]["mol2"])

    if "ssbond" not in setting["input"]["protein"] or setting["input"]["protein"]["ssbond"] is None:
        setting["input"]["protein"]["ssbond"] = []

    setting["general"]["yaml"] = YAML_PATH
    return setting
