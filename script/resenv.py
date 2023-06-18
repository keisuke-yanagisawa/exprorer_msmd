from typing import List, Set
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import gridData
from scipy .spatial import distance
from tqdm import tqdm
from joblib import Parallel, delayed
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
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
    interp = RegularGridInterpolator(dx.midpoints, dx.grid, method="nearest", fill_value=-1, bounds_error=False)
    values = interp(np.array(coords))

    if lt:
        resis = np.array(resis)[values < threshold]
    else:
        resis = np.array(resis)[values > threshold]

    return set(resis)


def __get_surrounded_resis_around_a_residue(model: Model,
                                            resi: int,
                                            env_distance: float
                                            ) -> Set[int]:
    residue_coords = uPDB.get_attr(model, "coord",
                                   sele=lambda a: uPDB.get_resi(a) == resi)
    environment_resis = set(
        uPDB.get_attr(model, "resid",
                      sele=lambda a: np.min(distance.cdist([a.get_coord()], residue_coords)) < env_distance)
    )
    return set(environment_resis)


def __wrapper(model_wo_water: Model,
              dx: gridData.Grid,
              focused_resname: str,
              res_atomnames: List[str],
              threshold: float,
              lt: bool,
              env_distance: float):
    focused_residue_resis = set(
        uPDB.get_attr(model_wo_water, "resid", sele=lambda a: uPDB.get_resname(a) == focused_resname)
    )

    resi_set = compute_SR_probe_resis(model_wo_water, dx, focused_resname, threshold, lt)

    ret_env_structs = []
    for resi in resi_set:
        exclude_resis = focused_residue_resis
        environment_resis = __get_surrounded_resis_around_a_residue(model_wo_water, resi, env_distance)
        environment_resis -= exclude_resis

        if len(environment_resis) == 0:
            continue

        sele = uPDB.Selector(lambda a: uPDB.get_resi(a) in (environment_resis | set([resi])))
        env_struct = uPDB.extract_substructure(model_wo_water, sele)
        ret_env_structs.append(env_struct)
    return ret_env_structs


def resenv(grid: gridData.Grid, ipdb: List[str], resn: str, res_atomnames: List[str], opdb: str,
           threshold: float = 0.2, lt: bool = False,
           env_distance: float = 4, verbose=False) -> Structure:
    """
    Extract probe which is on high-probability region with its environment (protein residues)
    """
    dx = grid

    out_helper = uPDB.PDBIOhelper(opdb)
    for path in ipdb:
        reader = uPDB.MultiModelPDBReader(path)

        # lst_of_lst = [__wrapper(model, dx, resn, res_atomnames, threshold, lt, env_distance)
        #               for model in tqdm(reader, desc="[extract res. env.]", disable=not verbose)]
        lst_of_lst = Parallel(n_jobs=1)(
            delayed(__wrapper)(model, dx, resn, res_atomnames, threshold, lt, env_distance)
            for model in tqdm(reader, desc="[extract res. env.]", disable=not verbose)
        )  # now parallel processing is disabled

        if lst_of_lst is None:
            continue

        for lst in lst_of_lst:
            for struct in lst:
                out_helper.save(struct)
    out_helper.close()  # it is necessary to read the output file below

    structs = uPDB.get_structure(opdb)
    if (len(structs) == 0):
        raise ValueError("No structures were extracted.")
    return structs
