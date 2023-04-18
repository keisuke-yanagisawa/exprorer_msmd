#!/usr/bin/python3

import os

from scipy import constants

from .utilities import util
from .utilities.executable import Cpptraj
from .utilities import GridUtil
from .utilities.Bio import PDB as uPDB
import gridData
import numpy as np

VERSION = "1.0.0"


def mask_generator(ref_struct, reference_grid, distance=None):
    mask = GridUtil.gen_distance_grid(reference_grid, ref_struct)
    # print(np.max(mask.grid), np.min(mask.grid), distance)
    if distance is not None:
        mask.grid = mask.grid < distance
    else:
        mask.grid = mask.grid < np.inf
    return mask


def convert_to_proba(g, mask_grid=None, normalize="snapshot", frames=1):
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


def convert_to_gfe(grid_path, mean_proba, temperature=300):
    pmap = gridData.Grid(grid_path)
    pmap.grid = np.where(pmap.grid <= 0, 1e-10, pmap.grid)  # avoid log(0)
    pmap.grid = -(constants.R / constants.calorie / constants.kilo) * temperature * np.log(pmap.grid / mean_proba)
    pmap.grid = np.where(pmap.grid > 3, 3, pmap.grid)  # Definition of GFE in the paper Raman et al., JCIM, 2013

    gfe_path = os.path.dirname(grid_path) + "/" \
        + "GFE" + "_" + os.path.basename(grid_path)
    pmap.export(gfe_path, type="double")

    pmap.grid = -pmap.grid
    invgfe_path = os.path.dirname(grid_path) + "/" \
        + "InvGFE" + "_" + os.path.basename(grid_path)
    pmap.export(invgfe_path, type="double")

    return gfe_path


def convert_to_pmap(grid_path, ref_struct, valid_distance, normalize="snapshot", frames=1):
    grid = gridData.Grid(grid_path)
    mask = mask_generator(ref_struct, grid, valid_distance)
    pmap = convert_to_proba(grid, mask.grid, frames=frames, normalize=normalize)

    pmap_path = os.path.dirname(grid_path) + "/" \
        + "PMAP" + "_" + os.path.basename(grid_path)
    pmap.export(pmap_path, type="double")
    return pmap_path


def parse_snapshot_setting(string):
    offset = "1"  # default parameter
    if len(string.split(":")) != 1:        # ofset is an option
        string, offset = string.split(":")
    start, stop = string.split("-")         # start and end is mandatory
    return start, stop, offset


def gen_pmap(dirpath, setting_general, setting_input, setting_pmap, traj, top, debug=False):

    traj_start, traj_stop, traj_offset \
        = parse_snapshot_setting(setting_pmap["snapshot"])

    name = setting_general["name"]

    trajectory = util.getabsolutepath(traj)
    topology = util.getabsolutepath(top)
    ref_struct = setting_input["protein"]["pdb"]
    probe_id = setting_input["probe"]["cid"]
    maps = setting_pmap["maps"]
    box_size = setting_pmap["map_size"]
    box_center = uPDB.get_structure(setting_input["protein"]["pdb"]).center_of_mass()

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
        maps=maps
    )

    pmap_paths = []
    for map in cpptraj_obj.maps:
        pmap_path = convert_to_pmap(map["grid"], ref_struct,
                                    setting_pmap["valid_dist"],
                                    frames=cpptraj_obj.frames,
                                    normalize=setting_pmap["normalization"])
        if setting_pmap["normalization"] == "GFE":
            struct_obj = uPDB.get_structure(ref_struct)
            protein_volume = uPDB.estimate_exclute_volume(struct_obj)
            mean_proba = map["num_probe_atoms"] / (cpptraj_obj.last_volume - protein_volume)
            pmap_path = convert_to_gfe(pmap_path, mean_proba,
                                       temperature=300)  # TODO: read temperature from setting
        pmap_paths.append(pmap_path)

    return pmap_paths
