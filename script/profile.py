import argparse
import tempfile
import os
from typing import List
import numpy as np
import numpy.typing as npt
import gridData

from script.utilities.Bio import PDB as uPDB
from script.utilities import const
from Bio import PDB
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


def __atomtype_select(struct: Structure,
                      atoms: List[str]):
    class AtomTypeSelect(PDB.Select):
        def __init__(self, atoms):
            self.atomtypes = atoms

        def accept_atom(self, atom):
            return uPDB.get_atom_attr(atom, "fullname") in self.atomtypes

    _, tmp = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)
    io = PDB.PDBIO()
    io.set_structure(struct)
    io.save(tmp, AtomTypeSelect(atoms))
    struct = uPDB.get_structure(tmp)
    os.remove(tmp)

    return struct


def get_all_atom_coordinates(struct: Structure) -> npt.NDArray[np.float64]:
    coords = []
    for m in struct:
        for a in m.get_atoms():
            coords.append(a.coord)
    return np.array(coords)


def get_all_atom_names(struct: Structure) -> npt.NDArray[np.str_]:
    names = []
    for m in struct:
        for a in m.get_atoms():
            names.append(uPDB.get_resname(a))
    return np.array(names)


def create_residue_interaction_profile(ipdb: str,
                                       atoms: List[str],
                                       residue_lst: List[str]) -> gridData.Grid:

    struct = uPDB.get_structure(ipdb)
    struct = __atomtype_select(struct, atoms)

    coords = get_all_atom_coordinates(struct)
    names = get_all_atom_names(struct)
    min_xyz = np.min(coords, axis=0)
    max_xyz = np.max(coords, axis=0)

    applicables = [s in residue_lst for s in names]
    target_coords = coords[np.where(applicables)]

    x_range = np.arange(np.floor(min_xyz[0]), np.ceil(max_xyz[0]) + 1, 1)
    y_range = np.arange(np.floor(min_xyz[1]), np.ceil(max_xyz[1]) + 1, 1)
    z_range = np.arange(np.floor(min_xyz[2]), np.ceil(max_xyz[2]) + 1, 1)

    hist, bins = np.histogramdd(target_coords, [x_range, y_range, z_range])

    return gridData.Grid(hist, edges=bins)

