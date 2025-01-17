#!/usr/bin/python3

from pathlib import Path
from typing import List, Union

import numpy as np
from gridData import Grid

from script.utilities.logger import logger

VERSION = "1.0.0"


def check(gs: List[Grid]) -> bool:
    """Check if all grids have same voxel sizes and positions.
    
    Args:
        gs: List of Grid objects to check
        
    Returns:
        True if all grids are compatible, False otherwise
    """
    return True  # TODO: Implement actual check


def grid_max(gs: List[Grid]) -> Grid:
    '''Generate maximum grid from list of grids.
    
    Args:
        gs: List of Grid objects to combine
        
    Returns:
        Grid: Grid containing maximum values across all input grids
        
    Raises:
        ValueError: If grids are incompatible
    '''
    if not check(gs):
        logger.error("ERROR: Grid(s) have different voxel sizes or/and positions")
        exit(1)

    ret = gs[0]
    ret.grid = np.max([g.grid for g in gs], axis=0)
    return ret


def gen_max_pmap(inpaths: List[Union[str, Path]], outpath: Union[str, Path]) -> str:
    '''Generate maximum probability map from multiple input maps.
    
    Args:
        inpaths: List of paths to input probability map files
        outpath: Path to save the output maximum probability map
        
    Returns:
        str: Path to the generated maximum probability map file
    '''
    gs = [Grid(n) for n in inpaths]
    max_pmap = grid_max(gs)
    max_pmap.export(outpath, type="double")
    return str(outpath)
