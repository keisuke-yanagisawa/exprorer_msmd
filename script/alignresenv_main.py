import argparse
from typing import List
from tqdm import tqdm
from Bio.PDB.Atom import Atom
from utilities.Bio import PDB as uPDB
from utilities.Bio.sklearn_interface import SuperImposer

DESCRIPTION = """
superimpose structures in accordance with specific atoms
"""


def align_res_env(ipdb: str,
                  resn: str,
                  opdb: str,
                  focused: List[str] = []):
    def selector(a: Atom):
        cond1 = uPDB.get_atom_attr(a, "resname") == resn
        cond2 = len(focused) == 0 or uPDB.get_atom_attr(a, "fullname") in focused
        return cond1 and cond2

    ref = uPDB.get_structure(ipdb[0])[0].copy()
    ref_probe_c_coords = uPDB.get_attr(ref, "coord", sele=selector)
    # print(ref_probe_c_coords)

    sup = SuperImposer()
    with uPDB.PDBIOhelper(opdb) as pdbio:
        for src in ipdb:
            struct = uPDB.get_structure(src)
            for model in tqdm(struct, desc="[align res. env.]", disable=not (VERBOSE or DEBUG)):

                # print(struct, i)
                probe_coords = uPDB.get_attr(model, "coord", sele=selector)
                sup.fit(probe_coords, ref_probe_c_coords)
                all_coords = uPDB.get_attr(model, "coord")
                uPDB.set_attr(model, "coord", sup.transform(all_coords))

                pdbio.save(model)
                # print(len(pdbio))
    return opdb


VERBOSE = None
DEBUG = None
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("-ipdb", required=True, nargs="+", help="input file path")
    parser.add_argument("-resn", required=True, help="residue name of the probe")
    parser.add_argument("-opdb", required=True, help="output file path")
    parser.add_argument("-v", dest="verbose", action="store_true")
    parser.add_argument("--focused", nargs="+", default=[],
                        help="please write ' CA ' or like that")

    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()
    VERBOSE = args.verbose
    DEBUG = args.debug
    align_res_env(args.ipdb, args.resn, args.opdb, args.focused)
