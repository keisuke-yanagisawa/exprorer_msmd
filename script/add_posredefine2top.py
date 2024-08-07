import os
from typing import List

import jinja2
import numpy.typing as npt
from scipy import constants

VERSION = "2.0.0"


def __position_restraint(atom_id_list: npt.ArrayLike, prefix: str, weight) -> str:
    """
    generate a string defining position restraint records
    """
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
    template = env.get_template("./template/position_restraints")
    return template.render(
        {
            "define_name": f"{prefix}{weight}",
            "weight": weight,
            "weight_in_calorie": weight * constants.calorie,
            "atom_id_list": atom_id_list,
        }
    )


def embed_posre(top_string: str, atom_id_list: npt.ArrayLike, prefix: str, strength: list[int]) -> str:
    """
    embed position restraint records into a given topology string
    """
    ret = []
    curr_section = None
    mol_count = 0
    for line in top_string.split("\n"):
        if line.startswith("["):
            curr_section = line[line.find("[") + 1 : line.find("]")].strip()
            if curr_section == "moleculetype":
                if mol_count == 1:
                    for s in strength:
                        ret.append(__position_restraint(atom_id_list, prefix, s))
                mol_count += 1

        ret.append(line)
    ret = "\n".join(ret)

    return ret
