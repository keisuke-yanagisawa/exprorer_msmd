#!/usr/bin/python3

import os
from pathlib import Path
from typing import Optional, Literal

import gridData
import numpy as np
import numpy.typing as npt
from scipy import constants

from script.utilities import GridUtil, util
from script.utilities.Bio import PDB as uPDB
from script.utilities.executable import Cpptraj

VERSION = "1.0.0"


def mask_generator(ref_struct: Path, reference_grid: gridData.Grid, distance: Optional[float] = None) -> gridData.Grid:
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
    g: gridData.Grid, mask_grid: Optional[npt.NDArray] = None, normalize: Literal["total", "snapshot"] = "snapshot", frames: int = 1
) -> gridData.Grid:
    if mask_grid is not None:
        values = g.grid[np.where(mask_grid)]
        # print(np.sum(g.grid), np.sum(values), np.where(mask_grid))
        if normalize == "snapshot":
            values /= frames
        elif normalize == "total":
            values /= np.sum(values)
        else:
            raise ValueError("Invalid normalization method")
        g.grid = np.full_like(g.grid, np.min([np.min(values), -1]))  # assign -1 for outside of mask
        g.grid[np.where(mask_grid)] = values
    else:
        g.grid /= np.sum(g.grid)
    return g


def convert_to_gfe(grid_path: str, mean_proba: float, temperature: float = 300) -> str:
    pmap = gridData.Grid(grid_path)
    pmap.grid = np.where(pmap.grid <= 0, 1e-10, pmap.grid)  # avoid log(0)
    pmap.grid = -(constants.R / constants.calorie / constants.kilo) * temperature * np.log(pmap.grid / mean_proba)
    pmap.grid = np.where(pmap.grid > 3, 3, pmap.grid)  # Definition of GFE in the paper Raman et al., JCIM, 2013

    gfe_path = os.path.dirname(grid_path) + "/" + "GFE" + "_" + os.path.basename(grid_path)
    pmap.export(gfe_path, type="double")

    pmap.grid = -pmap.grid
    invgfe_path = os.path.dirname(grid_path) + "/" + "InvGFE" + "_" + os.path.basename(grid_path)
    pmap.export(invgfe_path, type="double")

    return gfe_path


def convert_to_pmap(
    grid_path: Path, ref_struct: Path, valid_distance: float, normalize: Literal["total", "snapshot"] = "snapshot", frames: int = 1
):
    grid = gridData.Grid(grid_path)
    mask = mask_generator(ref_struct, grid, valid_distance)
    pmap = convert_to_proba(grid, mask.grid, frames=frames, normalize=normalize)

    pmap_path = os.path.dirname(grid_path) + "/" + "PMAP" + "_" + os.path.basename(grid_path)
    pmap.export(pmap_path, type="double")
    return pmap_path


def parse_snapshot_setting(string: str):
    if not string or ":" not in string and "-" not in string:
        raise ValueError("Invalid format. Expected 'start-stop' or 'start-stop:offset'")

    offset = "1"  # default parameter
    if ":" in string:  # offset is an option
        if string.count(":") > 1:
            raise ValueError("Too many ':' characters")
        string, offset = string.split(":")
        if not offset.isdigit():
            raise ValueError("Offset must be a positive integer")

    if "-" not in string:
        raise ValueError("Missing '-' separator")
    if string.count("-") > 1:
        raise ValueError("Too many '-' characters")
    
    try:
        start, stop = string.split("-")
    except ValueError:
        raise ValueError("Invalid format. Expected 'start-stop'")

    if not start or not stop:
        raise ValueError("Start and stop values must not be empty")
    if not start.isdigit() or not stop.isdigit():
        raise ValueError("Start and stop must be positive integers")
    if int(start) > int(stop):
        raise ValueError("Start frame must be less than or equal to stop frame")
    if int(offset) < 1:
        raise ValueError("Offset must be greater than 0")

    return start, stop, offset


def gen_pmap(
    dirpath: Path, setting_general: dict, setting_input: dict, setting_pmap: dict, traj: Path, top: Path, debug=False
):

    traj_start, traj_stop, traj_offset = parse_snapshot_setting(setting_pmap["snapshot"])

    name: str = setting_general["name"]

    trajectory = util.getabsolutepath(traj)
    topology = util.getabsolutepath(top)
    ref_struct = Path(setting_input["protein"]["pdb"])
    probe_id: str = setting_input["probe"]["cid"]
    maps: list = setting_pmap["maps"]
    box_size: int = setting_pmap["map_size"]
    box_center: npt.NDArray[np.float_] = uPDB.get_attr(
        uPDB.get_structure(setting_input["protein"]["pdb"]), "coord"
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
            normalize=setting_pmap["normalization"] if setting_pmap["normalization"] != "GFE" else "snapshot"
        )
        if setting_pmap["normalization"] == "GFE":
            struct_obj = uPDB.get_structure(ref_struct)
            protein_volume = uPDB.estimate_exclute_volume(struct_obj)
            mean_proba = map["num_probe_atoms"] / (cpptraj_obj.last_volume - protein_volume)
            pmap_path = convert_to_gfe(pmap_path, mean_proba, temperature=300)  # TODO: read temperature from setting
        pmap_paths.append(pmap_path)

    return pmap_paths
