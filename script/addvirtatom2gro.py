
import tempfile
from typing import List

from .utilities import gromacs
from .utilities.logger import logger
import numpy.typing as npt
import numpy as np
VERSION = "1.0.0"


def center_of_mass(atoms: List[gromacs.Gro_atom]) -> npt.NDArray:
    tot_weight = 0
    tot_coordinates = np.array([0, 0, 0])
    for a in atoms:
        tot_weight += a.atomic_mass
        tot_coordinates += a.atomic_mass * a.point
    return tot_coordinates / tot_weight


def addvirtatom2gro(gro_string: str,
                    probe_id: str):
    with tempfile.TemporaryDirectory() as tmpdirpath:
        with open(f"{tmpdirpath}/tmp.gro", "w") as fout:
            fout.write(gro_string)
        gro = gromacs.Gro(f"{tmpdirpath}/tmp.gro")

    atoms = [a.resi for a in gro.get_atoms(resn=probe_id)]
    resi_set = set(atoms)
    for resi in resi_set:
        mol = gro.get_atoms(resi=resi)
        com = center_of_mass(mol)
        virt_atom = gromacs.Gro_atom()
        virt_atom.point = com
        virt_atom.atomtype = "VIS"
        virt_atom.resn = mol[0].resn
        virt_atom.resi = resi
        gro.add_atom(virt_atom)

    logger.info(f"the system has {gro.molar(probe_id):.3f} M of {probe_id} cosolvents")
    return str(gro)

# import argparse
# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="add virtual atoms to a gro file")
#     parser.add_argument("-v,--version", action="version", version=VERSION)
#     parser.add_argument("-i", required=True, help="input gro file")
#     parser.add_argument("-o", required=True, help="output gro file")
#     parser.add_argument("-cname", required=True, help="cosolvent res name")
#     parser.add_argument("-vname", default="VIS", help="virtual interaction site atom name")
#     args = parser.parse_args()

#     infile = args.i
#     outfile = args.o

#     with open(infile) as fin:
#         gro_string = fin.read()
#     gro_string = addvirtatom2gro(gro_string, args.cname)
#     with open(outfile, "w") as fout:
#         fout.write(gro_string)
