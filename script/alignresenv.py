import argparse
import tempfile
from typing import List
from tqdm import tqdm
from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from script.utilities.Bio import PDB as uPDB
from script.utilities.Bio.sklearn_interface import SuperImposer

DESCRIPTION = """
superimpose structures in accordance with specific atoms
"""


def align_res_env(structs: List[Structure],
                  reference: Model,
                  resn: str,
                  focused: List[str] = []):

    def selector(a: Atom):
        cond1 = uPDB.get_atom_attr(a, "resname") == resn
        cond2 = len(focused) == 0 or uPDB.get_atom_attr(a, "fullname") in focused
        return cond1 and cond2

    ref_probe_c_coords = uPDB.get_attr(reference, "coord", sele=selector)
    # print(ref_probe_c_coords)

    sup = SuperImposer()
    _, tmppdb = tempfile.mkstemp(suffix=".pdb")
    with uPDB.PDBIOhelper(tmppdb) as pdbio:
        for struct in structs:
            for model in tqdm(struct, desc="[align res. env.]", disable=not (VERBOSE or DEBUG)):

                # print(struct, i)
                probe_coords = uPDB.get_attr(model, "coord", sele=selector)
                sup.fit(probe_coords, ref_probe_c_coords)
                all_coords = uPDB.get_attr(model, "coord")
                uPDB.set_attr(model, "coord", sup.transform(all_coords))

                pdbio.save(model)
                # print(len(pdbio))

    return uPDB.get_structure(tmppdb)


VERBOSE = None
DEBUG = None
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("-ipdb", required=True, nargs="+", help="input file path")
    parser.add_argument("-resn", required=True, help="residue name of the probe")
    parser.add_argument("-opdb", required=True, help="output file path")
    parser.add_argument("-v", dest="verbose", action="store_true")
    parser.add_argument("--focused", nargs="+", default=[],
                        help="focused atoms for probe alignment. please write ' CA ' or like that")

    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()
    VERBOSE = args.verbose
    DEBUG = args.debug

    ref_struct = uPDB.get_structure(args.ipdb[0])[0].copy()  # all structures are superimposed to this
    structures = [uPDB.get_structure(src) for src in args.ipdb]

    aligned = align_res_env(structures, ref_struct, args.resn, args.focused)
    with uPDB.PDBIOhelper(args.opdb) as pdbio:
        for model in aligned:
            pdbio.save(model)
