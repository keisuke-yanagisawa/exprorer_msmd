import argparse
from scipy import constants

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

def gen_atom_id_list(gro, target, RES, INV):
    atom_id_list = []

    MOLECULE = target == "molecule"
    
    mol_resi = -1
    mol_first_atom = -1 if MOLECULE else 1
    gro_strings = open(gro).readlines()[2:-2]
    for line in gro_strings:
        resi = int(line[:5])
        resn = line[5:10].strip()
        atom = line[10:15].strip()
        nr = int(line[15:20])
        if MOLECULE and (resn in RES) != INV:
            if mol_resi == -1:
                mol_resi = resi
                mol_first_atom = nr
            if mol_resi != resi:
                break
        if not atom.strip().startswith("H"):
            if (resn in RES) != INV:
                atom_id_list.append(nr-mol_first_atom+1)
    return atom_id_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="add position restraints for topology file")
    parser.add_argument("-v", action="store_true", help="invert selection")
    parser.add_argument("-res", required=True, nargs="+", help="residue names focused on")
    parser.add_argument("-target", choices=["protein", "molecule"], default="protein",
                        help="protein or molecule")
    parser.add_argument("-gro", required=True, help="input .gro file")
    parser.add_argument("-i", required=True, help="input topology file")
    parser.add_argument("-o", required=True, help="output topology file")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    atom_id_list = gen_atom_id_list(args.gro, args.target,
                                    args.res, args.v)
    with open(args.i) as fin:
        with open(args.o, "w") as fout:
            curr_section = None
            mol_count = 0
            for line in fin:
                if line.startswith("["):
                    curr_section = line[line.find("[")+1:line.find("]")].strip()
                    if curr_section == "moleculetype":
                        if mol_count == 1:
                            fout.write(position_restraint(atom_id_list))
                        mol_count += 1

                fout.write(line)
