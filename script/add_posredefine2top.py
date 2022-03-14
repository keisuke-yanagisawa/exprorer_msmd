import argparse
import os
import tempfile

from scipy import constants
import parmed as pmd

VERSION = "2.0.0"


def position_restraint(atom_id_list):
    ret_str = "; Position restraint\n"
    for weight in [1000, 500, 200, 100, 50, 20, 10, 0]:
        ret_str += f'#ifdef POSRES{weight}\n'
        ret_str += "[ position_restraints ]\n"
        ret_str += "; atom  type      fx      fy      fz\n"

        for atom_id in atom_id_list:
            c = int(constants.calorie * weight)
            ret_str += f"{atom_id: 6d}{1: 6d}{c: 6d}{c: 6d}{c: 6d}\n"
        ret_str += '#endif\n'
    return "\n"+ret_str+"\n"

def gen_protein_heavyatom_id_list(top_path):
    top = pmd.load_file(top_path)

    protein_resis = set(r.idx for r in top["@CA"].residues)
    atom_id_list = []
    for atom in top.atoms:
        if atom.element_name == "H":
            continue
        if atom.residue.idx not in protein_resis:
            continue
        atom_id_list.append(atom.idx+1)
    return atom_id_list

def add_posredefine2top(top_string):
    _, tmp = tempfile.mkstemp(suffix=".top")
    open(tmp, "w").write(top_string)
    atom_id_list = gen_protein_heavyatom_id_list(tmp)
    os.system(f"rm {tmp}")

    ret = []
    curr_section = None
    mol_count = 0
    for line in top_string.split("\n"):
        if line.startswith("["):
            curr_section = line[line.find("[")+1:line.find("]")].strip()
            if curr_section == "moleculetype":
                if mol_count == 1:
                    ret.append(position_restraint(atom_id_list))
                mol_count += 1
        ret.append(line+"\n")
    return "".join(ret)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="add position restraints for topology file")
    parser.add_argument("-i", required=True, help="input topology file")
    parser.add_argument("-o", required=True, help="output topology file")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    top_string = open(args.i).read()
    ret = add_posredefine2top(top_string)
    open(args.o, "w").write(ret)
    
