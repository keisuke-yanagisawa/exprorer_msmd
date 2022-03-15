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

def convert_to_pelmap(grid_path, ref_struct, frames, valid_distance):
    grid = gridData.Grid(grid_path)
    grid.grid /= frames
    mask = mask_generator(ref_struct, grid, valid_distance)
    grid.grid[~mask.grid] = -1
    
    pelmap_path = os.path.dirname(grid_path) + "/" \
                + "PELMAP" + "_" + os.path.basename(grid_path)
    grid.export(pelmap_path, type="double")
    return pelmap_path

def parse_snapshot_setting(string):
    offset = "1" # default parameter
    target, other = string.split("|")     # target name is mandatory
    if len(other.split(":")) != 1:        # ofset is an option
        other, offset = other.split(":")
    start, stop = other.split("-")         # start and end is mandatory
    return target, start, stop, offset

def gen_pelmap(index, setting_general, setting_input, setting_pelmap, debug=False):

    traj_target, traj_start, traj_stop, traj_offset \
        = parse_snapshot_setting(setting_pelmap["snapshot"])

    name = setting_general["name"]
    syspathdir = f"{setting_general['workdir']}/system{index}"

    trajectory = util.getabsolutepath(syspathdir) + f"/simulation/{traj_target}.xtc"
    topology   = util.getabsolutepath(syspathdir) + f"/top/{name}.top" # TODO: TEST_PROJECT should be "PREFIX"
    ref_struct = setting_input["protein"]["pdb"]
    probe_id   = setting_input["probe"]["cid"]
    maps       = setting_pelmap["maps"]
    box_size   = setting_pelmap["map_size"]

    cpptraj_obj = Cpptraj(debug=debug)
    cpptraj_obj.set(topology, trajectory, ref_struct, probe_id)
    cpptraj_obj.run(
        basedir=syspathdir, 
        prefix=name,
        box_size=box_size,
        traj_start=traj_start, 
        traj_stop=traj_stop, 
        traj_offset=traj_offset,
        maps = maps
    )
    
    pelmap_paths = []
    for grid_path in cpptraj_obj.grids:
        pelmap_path = convert_to_pelmap(grid_path, ref_struct, 
            cpptraj_obj.frames, setting_pelmap["valid_dist"])
        pelmap_paths.append(pelmap_path)
    
    return pelmap_paths