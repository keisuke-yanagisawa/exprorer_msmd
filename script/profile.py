from typing import List, Optional, Tuple, Union

import gridData
import numpy as np
import numpy.typing as npt
from Bio.PDB.Structure import Structure

from script.utilities.Bio import PDB as uPDB


class ProfileAnalyzer:
    """Analyzer for residue interaction profiles."""
    
    def __init__(self, structure: Structure) -> None:
        """Initialize analyzer with molecular structure.
        
        Args:
            structure: Structure containing aligned environments
        """
        self.structure = structure
        
    def _calc_minimum_bounding_box(self, coords: npt.NDArray[np.float64]) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Calculate minimum bounding box of coordinates.
        
        Args:
            coords: 3D coordinates array
            
        Returns:
            Tuple of minimum and maximum coordinates
        """
        min_xyz = np.min(coords, axis=0)
        max_xyz = np.max(coords, axis=0)
        return min_xyz, max_xyz
        
    def create_interaction_profile(
        self,
        target_residue_atoms: List[Tuple[str, str]]
    ) -> gridData.Grid:
        """Create residue interaction profile.
        
        Args:
            target_residue_atoms: List of (residue, atom) pairs to include
            
        Returns:
            Grid containing interaction profile
            
        Raises:
            ValueError: If no atoms found under specified names
        """
        sele = uPDB.Selector(
            lambda a: (uPDB.get_atom_attr(a, "resname"), uPDB.get_atom_attr(a, "fullname")) in target_residue_atoms
        )
        struct = uPDB.extract_substructure(self.structure, sele)
        if len(struct) == 0:
            raise ValueError("No atom found in the structure under the specified atom names")

        coords = uPDB.get_attr(struct, "coord")
        min_xyz, max_xyz = self._calc_minimum_bounding_box(coords)
        x_range = np.arange(np.floor(min_xyz[0]), np.ceil(max_xyz[0]) + 1, 1)
        y_range = np.arange(np.floor(min_xyz[1]), np.ceil(max_xyz[1]) + 1, 1)
        z_range = np.arange(np.floor(min_xyz[2]), np.ceil(max_xyz[2]) + 1, 1)

        hist, bins = np.histogramdd(coords, [x_range, y_range, z_range])
        return gridData.Grid(hist, edges=bins)


def create_residue_interaction_profile(struct: Structure, target_residue_atoms: List[Tuple[str, str]]) -> gridData.Grid:
    """Create residue interaction profile (wrapper for backward compatibility).
    
    Args:
        struct: Structure containing aligned environments
        target_residue_atoms: List of (residue, atom) pairs to include
            ex: [("ALA", " CA "), ("ALA", " CB "), ("ARG", " CB ")]
            
    Returns:
        Grid containing interaction profile
    """
    analyzer = ProfileAnalyzer(struct)
    return analyzer.create_interaction_profile(target_residue_atoms)
