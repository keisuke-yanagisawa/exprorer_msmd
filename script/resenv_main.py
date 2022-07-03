import numpy as np
import argparse
import io
import gridData
from scipy .spatial import distance
from tqdm import tqdm
from joblib import Parallel, delayed
import tempfile
from Bio import PDB
from utilities.Bio import PDB as uPDB

VERSION = "0.3.0"
DESCRIPTION = """
Extract probe which is on high-probability region 
with its environment
Output: residue environments
"""


def compute_SR_probe_resis(model, dx, resn, threshold, lt=False):
    """
    SR : Specified region
    return: set()
    """
    resis  = uPDB.get_attr(model, "resid",
                      sele=lambda a: uPDB.get_resname(a) == resn and not uPDB.is_hydrogen(a))
    coords = uPDB.get_attr(model, "coord",
                      sele=lambda a: uPDB.get_resname(a) == resn and not uPDB.is_hydrogen(a))
    XX, YY, ZZ = np.array(coords).T
    # print(coords, XX, YY, ZZ)
    values = dx.interpolated(XX, YY, ZZ)

    if lt:
        resis = np.array(resis)[values < threshold]
    else:
        resis = np.array(resis)[values > threshold]

    return set(resis)


class Selector(PDB.Select):
    def __init__(self, sele):
        self.sele = sele

    def accept_atom(self, atom):
        return self.sele(atom)

def wrapper(model, dx, resn, threshold, lt, env_distance):
    water_resis = set(
        uPDB.get_attr(model, "resid", sele=uPDB.is_water)
    )

    resi_set = compute_SR_probe_resis(model, dx, resn, threshold, lt)

    ret_env_structs = []
    for resi in resi_set:
        probe_coords = uPDB.get_attr(model, "coord",
                                      sele=lambda a: uPDB.get_resi(a) == resi)
        environment_resis = set(
            uPDB.get_attr(model, "resid",
                           sele=lambda a: np.min(distance.cdist([a.get_coord()], probe_coords)) < env_distance
                           and uPDB.get_resname(a) != resn)
        )

        if len(environment_resis - water_resis) == 0:
            continue

        pdbio = PDB.PDBIO()
        pdbio.set_structure(model)
        sele = Selector(lambda a: uPDB.get_resi(a) in (environment_resis|set([resi]) ) )

        with tempfile.NamedTemporaryFile(suffix=".pdb") as fp:
            pdbio.save(fp.name, select=sele)
            tmp = uPDB.get_structure(fp.name)[0]
            ret_env_structs.append(tmp)
    return ret_env_structs

def resenv(grid, ipdb, resn, opdb, 
           threshold=0.2, lt=False, env_distance=4, n_jobs=1):
    dx = gridData.Grid(grid)

    out_helper = uPDB.PDBIOhelper(opdb)
    for path in ipdb:
        reader = uPDB.MultiModelPDBReader(path)
        
        lst_of_lst = Parallel(n_jobs=n_jobs)(
            delayed(wrapper)(model, dx, resn, threshold, lt, env_distance)
            for model in tqdm(reader, desc="[extract res. env.]", disable=not (VERBOSE or DEBUG))
        )


        for lst in lst_of_lst:
            for struct in lst:
                out_helper.save(struct)

VERBOSE = None
DEBUG = None
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("-grid", required=True, help="input max-PMAP")
    parser.add_argument("-ipdb", required=True, nargs="+", help="input snapshot pdb file")
    parser.add_argument("-resn", required=True, help="probe residue name")
    parser.add_argument("-opdb", required=True, help="output pdb file")
    parser.add_argument("-lt,--less-than", dest="lt", action="store_true",
                        help="extract environment being less than threshold instead of greater than")
    parser.add_argument("-v", dest="verbose", action="store_true")

    parser.add_argument("--threshold", default=0.2, type=float,
                        help="probability threshold") # preferable protein surface
    parser.add_argument("--env-distance", default=4, type=int,
                        help="radius of the environment (distance from probe atoms)")
    parser.add_argument("--njobs", default=1, type=int)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    VERBOSE = args.verbose # assign to a global variable
    DEBUG   = args.debug   # assign to a global variable
    resenv(args.grid, args.ipdb, args.resn, args.opdb,
           args.threshold, args.lt, args.env_distance, args.njobs)