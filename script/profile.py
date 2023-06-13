import itertools
from typing import List, Tuple
import numpy as np
import gridData

from script.utilities.Bio import PDB as uPDB
from Bio.PDB.Structure import Structure

DESCRIPTION = """
generate residue interaction profile of a probe of interest
"""

residue_type_dict = {
    "anion": ["ASP", "GLU"],
    "cation": ["LYS", "HIE", "ARG"],
    "aromatic": ["TYR", "TRP", "PHE"],
    "hydrophilic": ["ASN", "GLN", "SER", "THR", "CYS"],
    "hydrophobic": ["ILE", "LEU", "VAL", "ALA", "PRO", "MET"],
    "neutral": ["GLY"],
    "gly": ["GLY"],
    "met": ["MET"]
}


def create_residue_interaction_profile(struct: Structure,
                                       target_residue_atoms: List[Tuple[str, str]]
                                       ) -> gridData.Grid:

    sele = uPDB.Selector(lambda a: (uPDB.get_atom_attr(a, "fullname"),
                                    uPDB.get_atom_attr(a, "resname")) in target_residue_atoms)
    struct = uPDB.extract_substructure(struct, sele)
    if len(struct) == 0:
        raise ValueError("No atom found in the structure under the specified atom names")

    coords = uPDB.get_attr(struct, "coord")
    min_xyz = np.min(coords, axis=0)
    max_xyz = np.max(coords, axis=0)
    x_range = np.arange(np.floor(min_xyz[0]), np.ceil(max_xyz[0]) + 1, 1)
    y_range = np.arange(np.floor(min_xyz[1]), np.ceil(max_xyz[1]) + 1, 1)
    z_range = np.arange(np.floor(min_xyz[2]), np.ceil(max_xyz[2]) + 1, 1)

    hist, bins = np.histogramdd(coords, [x_range, y_range, z_range])

    return gridData.Grid(hist, edges=bins)
