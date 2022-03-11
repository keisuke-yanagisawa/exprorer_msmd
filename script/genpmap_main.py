#!/usr/bin/python3

import argparse
import configparser
import tempfile
from os.path import basename, dirname
import os
from subprocess import getoutput as gop

from scipy import constants
import jinja2

from utilities.pmd import convert as pmd_convert
from utilities.util import expandpath
from utilities.executable import Cpptraj
from utilities import const, GridUtil
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
    
    pmap_path = grid_path + "_pmap.dx" # TODO: terrible naming
    pmap.export(pmap_path, type="double")
    return pmap_path


def gen_pmap(basedir, prot_param, cosolv_param, valid_distance, debug=False):
    params = configparser.ConfigParser()
    params.read(expandpath(prot_param), "UTF-8")
    params.read(expandpath(cosolv_param), "UTF-8")
    if "ReferenceStructure" not in params:
        params["ReferenceStructure"] = params["Protein"]

    trajectory = expandpath(basedir) + "/" + "simulation" + "/" + "TEST_PROJECT.xtc" # TODO: TEST_PROJECT should be "PREFIX"
    topology   = expandpath(basedir) + "/" + "top" + "/" + "TEST_PROJECT.top" # TODO: TEST_PROJECT should be "PREFIX"
    ref_struct = dirname(expandpath(prot_param))+"/"+basename(params["ReferenceStructure"]["pdb"])
    probe_id   = params["Cosolvent"]["cid"]

    cpptraj_obj = Cpptraj(debug=debug)
    cpptraj_obj.set(topology, trajectory, ref_struct, probe_id)
    cpptraj_obj.run(
        basedir=basedir, 
        prefix="test"
    )
    
    pmap_paths = []
    for grid_path in cpptraj_obj.grids:
        pmap_path = convert_to_pmap(grid_path, ref_struct, valid_distance)
        pmap_paths.append(pmap_path)
    
    return pmap_paths

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="generate PMAPs")
    parser.add_argument("-basedir", required=True,
                        help="objective directory")
    parser.add_argument("prot_param")
    parser.add_argument("cosolv_param")

    parser.add_argument("-d,--distance-threshold", dest="d", metavar="d", default=5, type=int,
                        help="distance from protein atoms.")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    pmap_paths = gen_pmap(args.basedir, args.prot_param, args.cosolv_param, args.d, args.debug)
