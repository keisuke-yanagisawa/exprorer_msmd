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

def convert_to_proba(g, mask_grid=None):
    if mask_grid is not None:
        values = g.grid[np.where(mask_grid)]
        # print(np.sum(g.grid), np.sum(values), np.where(mask_grid))
        values /= np.sum(values)
        g.grid = np.full_like(g.grid, np.min([np.min(values), -1])) #TODO: I could not remember why it is needed...??
        g.grid[np.where(mask_grid)] = values
    else:
        g.grid /= np.sum(g.grid)
    return g

def convert_to_pmap(grid_path, ref_struct, valid_distance):
    grid = gridData.Grid(grid_path)
    mask = mask_generator(ref_struct, grid, valid_distance)
    pmap = convert_to_proba(grid, mask.grid)
    
    pmap_path = os.path.dirname(grid_path) + "/" \
                + "PMAP" + "_" + os.path.basename(grid_path)
    pmap.export(pmap_path, type="double")
    return pmap_path


def gen_pmap(index, setting, debug=False):

    name = setting["general"]["name"]
    syspathdir = f"{setting['general']['workdir']}/system{index}"

    trajectory = util.getabsolutepath(syspathdir) + f"/simulation/{name}.xtc"
    topology   = util.getabsolutepath(syspathdir) + f"/prep/{name}.top" # TODO: TEST_PROJECT should be "PREFIX"
    ref_struct = setting["input"]["protein"]["pdb"]
    probe_id   = setting["input"]["probe"]["cid"]

    cpptraj_obj = Cpptraj(debug=debug)
    cpptraj_obj.set(topology, trajectory, ref_struct, probe_id)
    cpptraj_obj.run(
        basedir=syspathdir, 
        prefix=name
    )
    
    pmap_paths = []
    for grid_path in cpptraj_obj.grids:
        pmap_path = convert_to_pmap(grid_path, ref_struct, 
            setting["exprorer_msmd"]["pmap"]["valid_dist"])
        pmap_paths.append(pmap_path)
    
    return pmap_paths

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="generate PMAPs")
    parser.add_argument("-index", required=True, type=int)
    parser.add_argument("setting_yaml", help="yaml file for the MSMD")

    parser.add_argument("-v,--verbose", dest="verbose", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    if args.debug:
        logger.setLevel("debug")
    elif args.verbose:
        logger.setLevel("info")
    # else: logger level is "warn"

    setting = util.parse_yaml(args.setting_yaml)
    logger.info("PMAP generation")
    pmap_paths = gen_pmap(args.index, setting, args.debug)
