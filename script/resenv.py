from typing import List
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import gridData
from scipy .spatial import distance
from tqdm import tqdm
from joblib import Parallel, delayed
import tempfile
from Bio import PDB
from Bio.PDB.Model import Model
from script.utilities.Bio import PDB as uPDB

VERSION = "0.3.0"
DESCRIPTION = """
Extract probe which is on high-probability region
with its environment
Output: residue environments
"""


def compute_SR_probe_resis(model: Model,
                           dx: gridData.Grid,
                           resn: str,
                           threshold: float,
                           lt: bool = False):
    """
    This function enumerates the residue numbers `resis` of probe molecules 
    located in regions that exceed the `threshold` value.
    ------
    input:
        model: Model
            A snapshot of a molecular dynamics simulation
            containing probe molecules
        dx: gridData.Grid,
            A grid data containing values of a property of interest
        resn: str,
            A residue name of probe molecules
        threshold: float,
            A threshold value to determine the regions
        lt: bool = False
            If True, the regions with values less than the threshold are selected
    output:
        resis: set
            A set of residue numbers of probe molecules
    """
    resis = uPDB.get_attr(model, "resid",
                          sele=lambda a: uPDB.get_resname(a) == resn and not uPDB.is_hydrogen(a))
    coords = uPDB.get_attr(model, "coord",
                           sele=lambda a: uPDB.get_resname(a) == resn and not uPDB.is_hydrogen(a))
    # below interpolation is "3D-spline" interpolation. It shows undesirable behavior
    # XX, YY, ZZ = np.array(coords).T
    # values = dx.interpolated(XX, YY, ZZ)

    # Because of the reason, we changed to "nearest" interpolation to obtain stable results
    interp = RegularGridInterpolator(dx.midpoints, dx.grid, method="nearest", fill_value=-1, bounds_error=False)
    values = interp(np.array(coords))

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


def wrapper(model: Model,
            dx: gridData.Grid,
            resn: str,
            threshold: float,
            lt: bool,
            env_distance: float):
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
                          sele=lambda a: np.min(distance.cdist([a.get_coord()], probe_coords)) < env_distance and uPDB.get_resname(a) != resn)
        )

        if len(environment_resis - water_resis) == 0:
            continue

        pdbio = PDB.PDBIO()
        pdbio.set_structure(model)
        sele = Selector(lambda a: uPDB.get_resi(a) in (environment_resis | set([resi])))

        with tempfile.NamedTemporaryFile(suffix=".pdb") as fp:
            pdbio.save(fp.name, select=sele)
            tmp = uPDB.get_structure(fp.name)[0]
            ret_env_structs.append(tmp)
    return ret_env_structs


def resenv(grid: str, ipdb: List[str], resn: str, opdb: str,
           threshold: float = 0.2, lt: bool = False,
           env_distance: float = 4, n_jobs: int = 1, verbose=False):
    dx = gridData.Grid(grid)

    out_helper = uPDB.PDBIOhelper(opdb)
    for path in ipdb:
        reader = uPDB.MultiModelPDBReader(path)

        lst_of_lst = Parallel(n_jobs=n_jobs)(
            delayed(wrapper)(model, dx, resn, threshold, lt, env_distance)
            for model in tqdm(reader, desc="[extract res. env.]", disable=not verbose)
        )

        if lst_of_lst is None:
            continue

        for lst in lst_of_lst:
            for struct in lst:
                out_helper.save(struct)
