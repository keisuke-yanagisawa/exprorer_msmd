import argparse

from utilities import gromacs

VERSION = "1.0.0"


def center_of_mass(atoms):
    tot_weight = 0
    tot_coordinates = 0
    for a in atoms:
        tot_weight += a.atomic_mass
        tot_coordinates += a.atomic_mass * a.point
    return tot_coordinates/tot_weight


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="add virtual atoms to a gro file")
    parser.add_argument("-v,--version", action="version", version=VERSION)
    parser.add_argument("-i", required=True, help="input gro file")
    parser.add_argument("-o", required=True, help="output gro file")
    parser.add_argument("-cname", required=True, nargs="+", help="cosolvent res name")
    parser.add_argument("-vname", default="VIS", help="virtual interaction site atom name")
    args = parser.parse_args()

    infile = args.i
    outfile = args.o
    RES_NAMES = args.cname
    VIRT_NAME = args.vname

    gro = gromacs.Gro(infile)

    for RES_NAME in RES_NAMES:
        atoms = [a.resi for a in gro.get_atoms(resn=RES_NAME)]
        resi_set = set(atoms)
        for resi in resi_set:
            mol = gro.get_atoms(resi=resi)
            com = center_of_mass(mol)
            virt_atom = gromacs.Gro_atom()
            virt_atom.point = com
            virt_atom.atomtype = VIRT_NAME
            virt_atom.resn = mol[0].resn
            virt_atom.resi = resi
            gro.add_atom(virt_atom)
        print("the system has {:.3f} M of {} cosolvents"
              .format(gro.molar(RES_NAME), RES_NAME))

    fout = open(outfile, "w")
    fout.write("{}".format(gro))
    fout.close()
