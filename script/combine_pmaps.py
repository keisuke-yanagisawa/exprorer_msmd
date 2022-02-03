#! /usr/bin/python3

# standard library
import argparse
import os

# external libraries
import numpy as np
from scipy import stats
from gridData import Grid

# own library
from utilities import GridUtil


VERSION = "1.0.0"

def mask_generator(maskfile, referencefile=None, distance=5):
    if os.path.splitext(maskfile)[1] == ".pdb":
        mask = GridUtil.gen_distance_grid(referencefile, maskfile)
    else:
        mask = Grid(maskfile)

    mask.grid = mask.grid < distance
    return mask

def procedure_sum(gs, mask_arr=None):
    sum_g = np.sum(gs)
    if mask_arr is not None:
        values = sum_g.grid[np.where(mask_arr)]
        sum_g.grid = np.full_like(sum_g.grid, np.min([np.min(values), -1]))
        sum_g.grid[np.where(mask_arr)] = values
    return [sum_g]

def procedure_max(gs, mask_arr=None):
    ret = gs[0]
    ret.grid = np.max([g.grid for g in gs], axis=0)

    if mask_arr is not None:
        values = ret.grid[np.where(mask_arr)]
        ret.grid = np.full_like(ret.grid, np.min([np.min(values), -1]))
        ret.grid[np.where(mask_arr)] = values
    return [ret]

def procedure_proba(gs, mask_arr=None):
    ret = []
    for g in gs:
        if mask_arr is not None:
            values = g.grid[np.where(mask_arr)]
            values /= np.sum(values)
            g.grid = np.full_like(g.grid, np.min([np.min(values), -1]))
            g.grid[np.where(mask_arr)] = values
        else:
            g.grid /= np.sum(g.grid)
        ret.append(g)
    return ret

def procedure_zscore(gs, mask_arr=None):
    ret = []
    for g in gs:
        if mask_arr is not None:
            values = g.grid[np.where(mask_arr)]
            values = stats.zscore(values)
            g.grid = np.full_like(g.grid, np.min([np.min(values), -1]))
            g.grid[np.where(mask_arr)] = values
        else:
            g.grid = stats.zscore(g.grid, axis=None)
        ret.append(g)
    return ret

def procedure_none(gs, mask_arr=None):
    raise NotImplementedError("invalid mode was selected")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="combine multiple maps to a map")
    parser.add_argument("-i", metavar="IMAP", required=True, nargs="+",
                        help="input maps")
    parser.add_argument("-o", metavar="OMAP", required=True,
                        help="output map")
    parser.add_argument("--mode", default="sum",
                        choices=["sum", "max",
                                 "sum-probability", "max-probability", "probability-max",
                                 "sum-zscore", "max-zscore", "zscore-max"])
    parser.add_argument("-m", metavar="MASK",
                        help="mask dx file (storing distances to nearest protein atoms) or reference pdb file")
    parser.add_argument("--distance", metavar="d", default=5, type=int,
                        help="distance from protein atoms. it is valid only if mask is set")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    
    # generate mask: all true if mask is not provided
    tmp_grid = Grid(args.i[0])
    mask_arr = np.full_like(tmp_grid.grid, True, dtype=bool)
    if args.m is not None:
        mask_arr = mask_generator(args.m, args.i[0], args.distance).grid

    gs = [Grid(n) for n in args.i]
    sequence = args.mode.split("-")
    for proc_str in sequence:
        func_dict = {
            "sum": procedure_sum,
            "max": procedure_max,
            "probability": procedure_proba,
            "zscore": procedure_zscore
        }
        func = func_dict.get(proc_str,
                             procedure_none)
        gs = func(gs, mask_arr)

    assert len(gs) == 1
    ret = gs[0]

    # output processed map file
    ret.export(args.o, type="double")
