
import tempfile
from typing import List

from .utilities import gromacs
from .utilities.logger import logger
import numpy.typing as npt
import numpy as np
VERSION = "1.0.0"


def center_of_mass(atoms: List[gromacs.Gro_atom]) -> npt.NDArray:
    """
    Calculate the center of mass of a list of atoms
    """
    tot_weight = 0
    tot_coordinates = np.array([0.0, 0.0, 0.0])
    for a in atoms:
        tot_weight += a.atomic_mass
        tot_coordinates += a.atomic_mass * a.point
    return tot_coordinates / tot_weight


def addvirtatom2gro(gro_string: str,
                    probe_id: str) -> str:
    """
    Add virtual atoms to a gro file
    Virtual atoms are added to the center of mass of each probe
    """
    with tempfile.TemporaryDirectory() as tmpdirpath:
        with open(f"{tmpdirpath}/tmp.gro", "w") as fout:
            fout.write(gro_string)
        gro = gromacs.Gro(f"{tmpdirpath}/tmp.gro")

    atoms = [a.resi for a in gro.get_atoms(resn=probe_id)]
    resi_set = set(atoms)
    for resi in resi_set:
        mol = gro.get_atoms(resi=resi)
        com = center_of_mass(mol)
        virt_atom = gromacs.Gro_atom()
        virt_atom.point = com
        virt_atom.atomtype = "VIS"
        virt_atom.resn = mol[0].resn
        virt_atom.resi = resi
        gro.add_atom(virt_atom)

    logger.info(f"the system has {gro.molar(probe_id):.3f} M of {probe_id} cosolvents")
    return str(gro)
