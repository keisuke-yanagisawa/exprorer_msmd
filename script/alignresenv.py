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
                  focused: List[str] = [],
                  verbose: bool = False):

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
            for model in tqdm(struct, desc="[align res. env.]", disable=not verbose):

                # print(struct, i)
                probe_coords = uPDB.get_attr(model, "coord", sele=selector)
                sup.fit(probe_coords, ref_probe_c_coords)
                all_coords = uPDB.get_attr(model, "coord")
                uPDB.set_attr(model, "coord", sup.transform(all_coords))

                pdbio.save(model)
                # print(len(pdbio))

    return uPDB.get_structure(tmppdb)
