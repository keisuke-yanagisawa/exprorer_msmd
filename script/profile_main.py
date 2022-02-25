import argparse
import tempfile
import os
import numpy as np
from tqdm import tqdm
import gridData

from utilities.Bio import PDB as uPDB
from utilities import const
from Bio import PDB

DESCRIPTION = """
generate residue interaction profile of a probe of interest
"""

residue_type_dict = {
    "anion": ["ASP", "GLU"],
     "cation": ["LYS", "HIE", "ARG"],
     "aromatic": ["TYR", "TRP", "PHE"], 
     "hydrophilic": ["ASN", "GLN", "SER", "THR", "CYS"],
     "hydrophobic": ["ILE", "LEU", "VAL", "ALA", "PRO", "MET"],
     "neutral": ["GLY"],
     "gly": ["GLY"],
     "met": ["MET"]
}

def atomtype_select(struct, atoms):
    class AtomTypeSelect(PDB.Select):
        def __init__(self, atoms):
            self.atomtypes = atoms
        def accept_atom(self, atom):
            return uPDB.get_attr(atom, "fullname") in self.atomtypes

    _, tmp = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)
    io = PDB.PDBIO()
    io.set_structure(struct)
    io.save(tmp, AtomTypeSelect(atoms))
    struct = uPDB.get_structure(tmp)
    os.remove(tmp)

    return struct

def res_int_profile(ipdb, atoms, oprefix):
    struct = uPDB.get_structure(ipdb)
    struct = atomtype_select(struct, atoms)
    
    lim = np.array([[10000,-10000]]*3, dtype=float)
    coords = []
    names = []
    for m in tqdm(struct, desc="[generate res. int. profile]", disable=not (VERBOSE or DEBUG)):    
        for a in m.get_atoms():
            lim[:,0] = np.min([lim[:,0], a.coord], axis=0)
            lim[:,1] = np.max([lim[:,1], a.coord], axis=0)
            coords.append(a.coord)
            names.append(uPDB.get_resname(a))
    coords = np.array(coords)
    names  = np.array(names)
    
    for k, v in residue_type_dict.items():
        mesh_type = k
        applicables = [s in v for s in names]
        # print(np.array(applicables).shape)
        target_coords = coords[np.where(applicables)]
        # print(target_coords.shape)

        target_lim = np.zeros(lim.shape)
        target_lim[:,0] = np.floor(lim[:,0])
        target_lim[:,1] = np.ceil(lim[:,1]) + 1

        x_range = np.arange(*target_lim[0,:], 1)
        y_range = np.arange(*target_lim[1,:], 1)
        z_range = np.arange(*target_lim[2,:], 1)

        hist, bins = np.histogramdd(target_coords, [x_range, y_range, z_range])

        g = gridData.Grid(hist, edges=bins) 
        g.export(f"{oprefix}_{mesh_type}.dx", type="short")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("-ipdb", required=True, 
                        help="input file path")
    parser.add_argument("-oprefix", required=True, 
                        help="output file prefix")
    parser.add_argument("-atom", dest="atoms", 
                        help="atomtype to be focused on",
                        nargs="+", default=[" CB "])
    parser.add_argument("-v", dest="verbose", action="store_true")

    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()
    VERBOSE = args.verbose
    DEBUG = args.debug

    res_int_profile(args.ipdb, args.atoms, args.oprefix)