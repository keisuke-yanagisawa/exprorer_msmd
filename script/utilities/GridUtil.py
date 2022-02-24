import gridData
import numpy as np
import os
import errno
import copy
from scipy.spatial import distance
from Bio import PDB
from tqdm import tqdm

def gen_distance_grid(g_ref, pdbpath, verbose=True):
    grid_points = np.array([p for p in g_ref.centers()], dtype="int8")
    pdb = PDB.PDBParser().get_structure(pdbpath, pdbpath)
    coords = []
    for atom in pdb.get_atoms():
        if(atom.full_id[3][0].strip() != ""):
            continue
        coords.append(atom.coord)
    coords = np.array(coords, dtype="float32")
    
    # We change the implementation because of the memory usage
    # min_dist = distance.cdist(grid_points, coords).min(axis=1)
    min_dist = np.array([np.inf]*len(grid_points))
    for coord in tqdm(coords, desc="[gen_distance_grid]", disable=not verbose):
        min_dist = np.min([min_dist, np.linalg.norm(grid_points-coord, axis=1)], axis=0)


    distance_grid = copy.deepcopy(g_ref)
    distance_grid.grid = min_dist.reshape(distance_grid.grid.shape)
    return distance_grid
