from typing import List
from scipy import constants

VERSION = "2.0.0"


def position_restraint(atom_id_list: List[int], weight) -> str:
    """
    generate a string defining position restraint records
    """
    ret_str = "; Position restraint\n"
    ret_str += f'#ifdef POSRES{weight}\n'
    ret_str += "[ position_restraints ]\n"
    ret_str += "; atom  type      fx      fy      fz\n"

    for atom_id in atom_id_list:
        c = int(constants.calorie * weight)
        ret_str += f"{atom_id: 6d}{1: 6d}{c: 6d}{c: 6d}{c: 6d}\n"
    ret_str += '#endif\n'
    return ret_str


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
                atom_id_list: List[int]) -> str:
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
                    ret.append(position_restraint(atom_id_list, 1000))
                    ret.append(position_restraint(atom_id_list, 500))
                    ret.append(position_restraint(atom_id_list, 200))
                    ret.append(position_restraint(atom_id_list, 100))
                    ret.append(position_restraint(atom_id_list, 50))
                    ret.append(position_restraint(atom_id_list, 20))
                    ret.append(position_restraint(atom_id_list, 10))
                    ret.append(position_restraint(atom_id_list, 0))
                mol_count += 1

        ret.append(line)
    ret = "\n".join(ret)

    return ret


def add_posredefine2top(top_string: str,
                        gro_string: str,
                        cid: str) -> str:
    """
    add position restraint records to a given topology string
    """

    atom_id_list = gen_atom_id_list(
        gro_string,
        ["WAT", "Na+", "Cl-", "CA", "MG", "ZN", "CU", cid]
    )
    return embed_posre(top_string, atom_id_list)
