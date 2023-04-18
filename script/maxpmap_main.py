#!/usr/bin/python3

import argparse

from utilities.logger import logger
from gridData import Grid
import numpy as np

VERSION = "1.0.0"


def check(gs):
    return True


def grid_max(gs):
    # TODO: check all grids have the same voxel fields
    if not check(gs):
        logger.error("ERROR: Grid(s) have different voxel fields")
        exit(1)

    ret = gs[0]
    ret.grid = np.max([g.grid for g in gs], axis=0)
    return ret


def gen_max_pmap(inpaths, outpath):
    gs = [Grid(n) for n in inpaths]
    max_pmap = grid_max(gs)
    max_pmap.export(outpath, type="double")
    return outpath


DEBUG = None  # global variable
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="generate PMAPs")
    parser.add_argument("pmap_paths", nargs="+",
                        help="PMAPs will be integrated")
    parser.add_argument("mappmap_path",
                        help="destination path to output max-PMAP")

    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()
    DEBUG = args.debug  # assign to a global variable

    max_pmap_path = gen_max_pmap(args.pmap_paths, args.mappmap_path)
