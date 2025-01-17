#!/usr/bin/python3

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, TypedDict, Union, cast

import gridData
import numpy.typing as npt
from script.generate_msmd_system import InputSettings

class MapSettings(TypedDict):
    suffix: str
    selector: str

class PmapSettings(TypedDict):
    snapshot: str
    valid_dist: float
    map_size: int
    normalization: str
    maps: List[MapSettings]

class GeneralSettings(TypedDict):
    name: str
import numpy as np
import numpy.typing as npt
from scipy import constants

from script.utilities import GridUtil, util
from script.utilities.Bio import PDB as uPDB
from script.utilities.executable import Cpptraj

VERSION = "1.0.0"


def mask_generator(
    ref_struct: Path,
    reference_grid: gridData.Grid,
    distance: Optional[float] = None
) -> gridData.Grid:
    """
    input
        ref_struct: path to reference structure
        reference_grid: gridData.Grid object
        distance: distance threshold for mask
    output:
        mask: gridData.Grid object containing boolean values
    """
    mask = GridUtil.gen_distance_grid(reference_grid, ref_struct)
    # print(np.max(mask.grid), np.min(mask.grid), distance)
    if distance is not None:
        mask.grid = mask.grid < distance
    else:
        mask.grid = mask.grid < np.inf
    return mask


def convert_to_proba(
    g: gridData.Grid,
    mask_grid: Optional[npt.NDArray[np.float64]] = None,
    normalize: str = "snapshot",
    frames: int = 1
) -> gridData.Grid:
    if mask_grid is not None:
        values = g.grid[np.where(mask_grid)]
        # print(np.sum(g.grid), np.sum(values), np.where(mask_grid))
        if normalize == "snapshot" or normalize == "GFE":
            values /= frames
        else:
            values /= np.sum(values)
        g.grid = np.full_like(g.grid, np.min([np.min(values), -1]))  # assign -1 for outside of mask
        g.grid[np.where(mask_grid)] = values
    else:
        g.grid /= np.sum(g.grid)
    return g


class ProbabilityMap:
    """Class for generating and manipulating probability maps."""
    
    def __init__(self, grid: gridData.Grid) -> None:
        """Initialize probability map with grid data.
        
        Args:
            grid: Grid data object containing probability values
        """
        self.grid = grid
        self._mask: Optional[npt.NDArray[np.float64]] = None
        
    @classmethod
    def from_file(cls, path: Union[str, Path]) -> 'ProbabilityMap':
        """Create ProbabilityMap from a grid file.
        
        Args:
            path: Path to grid file
            
        Returns:
            New ProbabilityMap instance
        """
        return cls(gridData.Grid(str(path)))
        
    def generate_mask(self, ref_struct: Path, distance: Optional[float] = None) -> None:
        """Generate mask from reference structure.
        
        Args:
            ref_struct: Path to reference structure
            distance: Optional distance threshold for mask
        """
        mask_grid = mask_generator(ref_struct, self.grid, distance)
        self._mask = mask_grid.grid
        
    def convert_to_probability(self, frames: int = 1, normalize: str = "snapshot") -> None:
        """Convert grid values to probabilities.
        
        Args:
            frames: Number of frames used
            normalize: Normalization method ("snapshot", "GFE", or "total")
        """
        converted = convert_to_proba(self.grid, self._mask, normalize=normalize, frames=frames)
        self.grid = converted
        
    def save(self, path: Union[str, Path], type: str = "double") -> None:
        """Save probability map to file.
        
        Args:
            path: Output file path
            type: Data type for output
        """
        self.grid.export(str(path), type=type)
        
    def to_gfe(self, mean_proba: float, temperature: float = 300) -> None:
        """Convert probability map to Gibbs free energy.
        
        Args:
            mean_proba: Mean probability for normalization
            temperature: Temperature in Kelvin
        """
        self.grid.grid = np.where(self.grid.grid <= 0, 1e-10, self.grid.grid)
        self.grid.grid = -(constants.R / constants.calorie / constants.kilo) * temperature * np.log(self.grid.grid / mean_proba)
        self.grid.grid = np.where(self.grid.grid > 3, 3, self.grid.grid)
        
    def to_inverse_gfe(self) -> None:
        """Convert grid to inverse Gibbs free energy."""
        self.grid.grid = -self.grid.grid

def convert_to_gfe(
    grid_path: str,
    mean_proba: float,
    temperature: float = 300
) -> str:
    """Convert probability map to GFE (wrapper for backward compatibility).
    
    Args:
        grid_path: Path to probability map file
        mean_proba: Mean probability for normalization
        temperature: Temperature in Kelvin
        
    Returns:
        Path to generated GFE file
    """
    pmap = ProbabilityMap.from_file(grid_path)
    pmap.to_gfe(mean_proba, temperature)
    
    gfe_path = os.path.dirname(grid_path) + "/" + "GFE" + "_" + os.path.basename(grid_path)
    pmap.save(gfe_path)
    
    pmap.to_inverse_gfe()
    invgfe_path = os.path.dirname(grid_path) + "/" + "InvGFE" + "_" + os.path.basename(grid_path)
    pmap.save(invgfe_path)
    
    return gfe_path


def convert_to_pmap(
    grid_path: Path,
    ref_struct: Path,
    valid_distance: float,
    normalize: str = "snapshot",
    frames: int = 1
) -> str:
    """Convert grid to probability map (wrapper for backward compatibility).
    
    Args:
        grid_path: Path to input grid file
        ref_struct: Path to reference structure
        valid_distance: Distance threshold for mask
        normalize: Normalization method
        frames: Number of frames used
        
    Returns:
        Path to generated probability map file
    """
    pmap = ProbabilityMap.from_file(grid_path)
    pmap.generate_mask(ref_struct, valid_distance)
    pmap.convert_to_probability(frames=frames, normalize=normalize)
    
    pmap_path = os.path.dirname(grid_path) + "/" + "PMAP" + "_" + os.path.basename(grid_path)
    pmap.save(pmap_path)
    return pmap_path


def parse_snapshot_setting(string: str) -> Tuple[str, str, str]:
    offset = "1"  # default parameter
    if len(string.split(":")) != 1:  # ofset is an option
        string, offset = string.split(":")
    start, stop = string.split("-")  # start and end is mandatory
    return start, stop, offset


def gen_pmap(
    dirpath: Path,
    setting_general: GeneralSettings,
    setting_input: InputSettings,
    setting_pmap: PmapSettings,
    traj: Path,
    top: Path,
    debug: bool = False
) -> List[str]:

    traj_start, traj_stop, traj_offset = parse_snapshot_setting(setting_pmap["snapshot"])

    name: str = setting_general["name"]

    trajectory = util.getabsolutepath(traj)
    topology = util.getabsolutepath(top)
    ref_struct = Path(str(setting_input["protein"]["pdb"]))
    probe_id: str = cast(str, setting_input["probe"]["cid"])
    maps: List[MapSettings] = setting_pmap["maps"]
    box_size: int = cast(int, setting_pmap["map_size"])
    pdb_path = Path(str(setting_input["protein"]["pdb"]))
    box_center: npt.NDArray[np.float64] = uPDB.get_attr(
        uPDB.get_structure(pdb_path), "coord"
    ).mean(axis=0)
    # structure.center_of_mass() may return "[ nan nan nan ]" due to unspecified atomic weight

    cpptraj_obj = Cpptraj(debug=debug)
    cpptraj_obj.set(topology, trajectory, ref_struct, probe_id)
    cpptraj_obj.run(
        basedir=dirpath,
        prefix=name,
        box_size=box_size,
        box_center=box_center,
        traj_start=traj_start,
        traj_stop=traj_stop,
        traj_offset=traj_offset,
        maps=maps,
    )

    pmap_paths = []
    for map in cpptraj_obj.maps:
        pmap_path = convert_to_pmap(
            map["grid"],
            ref_struct,
            setting_pmap["valid_dist"],
            frames=cpptraj_obj.frames,
            normalize=setting_pmap["normalization"],
        )
        if setting_pmap["normalization"] == "GFE":
            struct_obj = uPDB.get_structure(ref_struct)
            protein_volume = uPDB.estimate_exclute_volume(struct_obj)
            mean_proba = map["num_probe_atoms"] / (cpptraj_obj.last_volume - protein_volume)
            pmap_path = convert_to_gfe(pmap_path, mean_proba, temperature=300)  # TODO: read temperature from setting
        pmap_paths.append(pmap_path)

    return pmap_paths
