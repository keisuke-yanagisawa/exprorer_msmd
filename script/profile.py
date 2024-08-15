from typing import List, Tuple

import gridData
import numpy as np
from Bio.PDB.Structure import Structure

from script.utilities.Bio import PDB as uPDB


def __calc_minimum_bounding_box(coords) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the minimum bounding box of a set of coordinates (3D)
    :param coords: 3D coordinates
    :return: min and max coordinates of the bounding box
    """
    min_xyz = np.min(coords, axis=0)
    max_xyz = np.max(coords, axis=0)
    return min_xyz, max_xyz


def create_residue_interaction_profile(struct: Structure, target_residue_atoms: tuple[tuple[str, str]]) -> gridData.Grid:
    """
    struct: Bio.PDB.Structure
        A structure containing multiple models of aligned environments
    target_residue_atoms: tuple[tuple[str, str]]
        A list of residue-atom pairs which are to be included in the profile
        ex: (("ALA", " CA "), ("ALA", " CB "), ("ARG", " CB "))
    """

    sele = uPDB.Selector(
        lambda a: (uPDB.get_atom_attr(a, "resname"), uPDB.get_atom_attr(a, "fullname")) in target_residue_atoms
    )
    struct = uPDB.extract_substructure(struct, sele)
    if len(struct) == 0:
        raise ValueError("No atom found in the structure under the specified atom names")

    coords = uPDB.get_attr(struct, "coord")
    min_xyz, max_xyz = __calc_minimum_bounding_box(coords)
    x_range = np.arange(np.floor(min_xyz[0]), np.ceil(max_xyz[0]) + 1, 1)
    y_range = np.arange(np.floor(min_xyz[1]), np.ceil(max_xyz[1]) + 1, 1)
    z_range = np.arange(np.floor(min_xyz[2]), np.ceil(max_xyz[2]) + 1, 1)

    hist, bins = np.histogramdd(coords, [x_range, y_range, z_range])

    return gridData.Grid(hist, edges=bins)
