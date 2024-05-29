import os
from typing import List
import jinja2
from scipy import constants
import numpy.typing as npt

VERSION = "2.0.0"


def position_restraint(atom_id_list: npt.ArrayLike, prefix: str, weight) -> str:
    """
    generate a string defining position restraint records
    """
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
    template = env.get_template("./template/position_restraints")
    return template.render({
        "define_name": f"{prefix}{weight}",
        "weight": weight,
        "weight_in_calorie": weight * constants.calorie,
        "atom_id_list": atom_id_list
    })


def gen_atom_id_list(gro_string: str,
                     except_resns: List[str]) -> List[int]:
    """
    generate a list of atom ids for a given residue name
    """
    atom_id_list = []
    mol_first_atom = 1
    gro_strings = gro_string.split("\n")[2:-2]
    for line in gro_strings:
        resn = line[5:10].strip()
        atom = line[10:15].strip()
        nr = int(line[15:20])
        if not atom.strip().startswith("H"):
            if not (resn in except_resns):
                atom_id_list.append(nr - mol_first_atom + 1)
    return atom_id_list


def embed_posre(top_string: str,
                atom_id_list: npt.ArrayLike,
                prefix: str,
                strength: list[int]) -> str:
    """
    embed position restraint records into a given topology string
    """
    ret = []
    curr_section = None
    mol_count = 0
    for line in top_string.split("\n"):
        if line.startswith("["):
            curr_section = line[line.find("[") + 1:line.find("]")].strip()
            if curr_section == "moleculetype":
                if mol_count == 1:
                    for s in strength:
                        ret.append(position_restraint(atom_id_list, prefix, s))
                mol_count += 1

        ret.append(line)
    ret = "\n".join(ret)

    return ret


def add_posredefine2top(top_string: str,
                        atom_id_list: list[int]) -> str:

    return embed_posre(top_string, atom_id_list)
