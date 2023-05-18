#!/usr/bin/python3

from typing import List

from script.utilities.logger import logger
from gridData import Grid
import numpy as np

VERSION = "1.0.0"


def check(gs: List[Grid]):
    return True


def grid_max(gs: List[Grid]):
    # TODO: check all grids have the same voxel fields
    if not check(gs):
        logger.error("ERROR: Grid(s) have different voxel fields")
        exit(1)

    ret = gs[0]
    ret.grid = np.max([g.grid for g in gs], axis=0)
    return ret


def gen_max_pmap(inpaths: List[str],
                 outpath: str):
    """
    Generate max pmap file from multiple pmap files
    Note that input/output is dx FILEs not grid objects
    """
    gs = [Grid(n) for n in inpaths]
    max_pmap = grid_max(gs)
    max_pmap.export(outpath, type="double")
    return outpath