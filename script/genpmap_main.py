#!/usr/bin/python3

import argparse
import configparser
import tempfile
from os.path import basename, dirname
import os
from subprocess import getoutput as gop

from scipy import constants
import jinja2

from .utilities.pmd import convert as pmd_convert
from .utilities.util import expandpath
from .utilities import util
from .utilities.executable import Cpptraj
from .utilities import const, GridUtil
from .utilities.logger import logger
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
        if normalize == "snapshot":
            values /= frames
        else:
            values /= np.sum(values)
        g.grid = np.full_like(g.grid, np.min([np.min(values), -1])) # assign -1 for outside of mask
        g.grid[np.where(mask_grid)] = values
    else:
        g.grid /= np.sum(g.grid)
    return g

def convert_to_pmap(grid_path, ref_struct, valid_distance, normalize="snapshot", frames=1):
    grid = gridData.Grid(grid_path)
    mask = mask_generator(ref_struct, grid, valid_distance)
    pmap = convert_to_proba(grid, mask.grid, frames=frames, normalize=normalize)
    
    pmap_path = os.path.dirname(grid_path) + "/" \
                + "PMAP" + "_" + os.path.basename(grid_path)
    pmap.export(pmap_path, type="double")
    return pmap_path

def parse_snapshot_setting(string):
    offset = "1" # default parameter
    if len(string.split(":")) != 1:        # ofset is an option
        string, offset = string.split(":")
    start, stop = string.split("-")         # start and end is mandatory
    return start, stop, offset

def gen_pmap(dirpath, setting_general, setting_input, setting_pmap, traj, top, debug=False):

    traj_start, traj_stop, traj_offset \
        = parse_snapshot_setting(setting_pmap["snapshot"])

    name = setting_general["name"]

    trajectory = util.getabsolutepath(traj)
    topology   = util.getabsolutepath(top)
    ref_struct = setting_input["protein"]["pdb"]
    probe_id   = setting_input["probe"]["cid"]
    maps       = setting_pmap["maps"]
    box_size  = setting_pmap["map_size"]

    cpptraj_obj = Cpptraj(debug=debug)
    cpptraj_obj.set(topology, trajectory, ref_struct, probe_id)
    cpptraj_obj.run(
        basedir=dirpath, 
        prefix=name,
        box_size=box_size,
        traj_start=traj_start, 
        traj_stop=traj_stop, 
        traj_offset=traj_offset,
        maps = maps
    )
    
    pmap_paths = []
    for grid_path in cpptraj_obj.grids:
        pmap_path = convert_to_pmap(grid_path, ref_struct, 
            setting_pmap["valid_dist"], 
            frames=cpptraj_obj.frames,
            normalize=setting_pmap["normalization"])
        pmap_paths.append(pmap_path)
    
    return pmap_paths