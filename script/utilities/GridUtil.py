import gridData
import numpy as np
import os
import errno
import copy
from scipy.spatial import distance
from Bio import PDB


def gen_distance_grid(g_ref, pdbpath):
    grid_points = np.array([p for p in g_ref.centers()], dtype="int8")
    pdb = PDB.PDBParser().get_structure(pdbpath, pdbpath)
    coords = []
    for atom in pdb.get_atoms():
        if(atom.full_id[3][0].strip() != ""):
            continue
        coords.append(atom.coord)
    coords = np.array(coords, dtype="float32")
    min_dist = distance.cdist(grid_points, coords).min(axis=1)


    distance_grid = copy.deepcopy(g_ref)
    distance_grid.grid = min_dist.reshape(distance_grid.grid.shape)
    return distance_grid
