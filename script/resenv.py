from typing import List, Optional, Set, Union
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import gridData
from scipy.spatial import distance
from tqdm import tqdm
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from script.utilities.Bio import PDB as uPDB

VERSION = "0.3.0"
DESCRIPTION = """
Extract probe which is on high-probability region
with its environment
Output: residue environments
"""


def compute_SR_probe_resis(model: Union[Structure, Model],
                           dx: gridData.Grid,
                           resn: str,
                           threshold: float,
                           lt: bool = False):
    """
    This function enumerates the residue numbers `resis` of probe molecules
    located in regions that exceed the `threshold` value.
    ------
    input:
        model: Union[Structure, Model]
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


def __get_surrounded_resis_around_a_residue(model: Union[Structure, Model],
                                            focused_resi: int,
                                            env_distance: float
                                            ) -> Set[int]:
    focused_residue_coords = uPDB.get_attr(model, "coord",
                                           sele=lambda a: uPDB.get_resi(a) == focused_resi)

    all_resis = uPDB.get_attr(model, "resid")
    all_coords = uPDB.get_attr(model, "coord")
    is_near_atom = np.min(distance.cdist(focused_residue_coords, all_coords), axis=0) < env_distance
    environment_resis = set(all_resis[is_near_atom])
    return set(environment_resis)


def __wrapper(model_wo_water: Union[Structure, Model],
              dx: gridData.Grid,
              focused_resname: str,
              res_atomnames: List[str] = [" CB "],
              threshold: float = 0.2,
              lt: bool = False,
              env_distance: float = 4.0) -> Optional[Structure]:
    focused_residue_resis = set(
        uPDB.get_attr(model_wo_water, "resid", sele=lambda a: uPDB.get_resname(a) == focused_resname)
    )
    # TODO: remove un-focusing atoms (not res_atomnames atoms)

    resi_set = compute_SR_probe_resis(model_wo_water, dx, focused_resname, threshold, lt)

    ret_env_structs = []
    for resi in resi_set:
        environment_resis = __get_surrounded_resis_around_a_residue(model_wo_water, resi, env_distance)
        environment_resis -= focused_residue_resis

        if len(environment_resis) == 0:
            # there is no protein residue around the probe molecule
            continue

        sele = uPDB.Selector(lambda a: uPDB.get_resi(a) in (environment_resis | set([resi])))
        env_struct = uPDB.extract_substructure(model_wo_water, sele)
        ret_env_structs.append(env_struct)

    ret = None
    if len(ret_env_structs) != 0:
        ret = uPDB.concatenate_structures(ret_env_structs)
    return ret


def resenv(grid: gridData.Grid,
           trajectory: uPDB.MultiModelPDBReader,
           resn: str,
           res_atomnames: List[str],
           threshold: float = 0.2,
           lt: bool = False,
           env_distance: float = 4,
           verbose: bool = False) -> Structure:
    """
    Extract probe which is on high-probability region with its environment (protein residues)
    """

    ret = []
    environments = [__wrapper(model, grid, resn, res_atomnames, threshold, lt, env_distance)
                    for model in tqdm(trajectory, desc="[extract res. env.]", disable=not verbose)]
    environments = [e for e in environments if e is not None]
    ret.extend(environments)

    if (len(ret) == 0):
        raise ValueError("No structures were extracted.")
    return uPDB.concatenate_structures(ret)
