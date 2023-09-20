import numpy as np
import numpy.typing as npt
from msmd.probe.parameter import Probe
from msmd.system import System
from msmd.standard_library.logging.logger import logger

import tempfile
from typing import List

from .simulation.gromacs.gro import Gro, GroAtom

VERSION = "1.0.0"


VERSION = "1.0.0"

VIS_INFO = """
[ atomtypes ]
VIS      VIS          0.00000  0.00000   V     0.00000e+00   0.00000e+00 ; virtual interaction site

[ nonbond_params ]
; i j func sigma epsilon
VIS   VIS    1  {sigma:1.6e}   {epsilon:1.6e}
"""


def addvirtatom2top(top_string: str,
                    probe_names: str,
                    sigma: float = 2,
                    epsilon: float = 4.184e-6) -> str:
    """
    Add definition of a virtual atom to a top file
    Pseudo repulsion term (VIS-VIS nonbond LJ parameter) is added to the top file.
    """

    ret = []
    curr_section = None
    now_mol = None
    natoms = 0
    for line in top_string.split("\n"):
        line = line.split(";")[0].strip()
        if line.startswith("["):
            prev_section = curr_section
            if prev_section == "atomtypes":
                ret.append(
                    VIS_INFO.format(sigma=sigma, epsilon=epsilon)
                )
            elif prev_section == "atoms":
                if now_mol == probe_names:  # TODO
                    ret.append(f"""
                    {natoms+1: 5d}        VIS      1    {now_mol}    VIS  {natoms+1: 5d} 0.00000000   0.000000
                    [ virtual_sitesn ]
                    {natoms+1: 5d}   2  {' '.join([str(x) for x in range(1, natoms+1)])}
                    """)
                natoms = 0
            curr_section = line[line.find("[") + 1:line.find("]")].strip()
            if curr_section == "moleculetype":
                now_mol = None
        elif curr_section == "atoms" and line != "":
            natoms += 1
        elif curr_section == "moleculetype" and now_mol is None and line != "":
            now_mol = line.split()[0].strip()

        ret.append(line)
    ret = "\n".join(ret)
    return ret


def center_of_mass(atoms: List[GroAtom]) -> npt.NDArray:
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
        gro = Gro(f"{tmpdirpath}/tmp.gro")

    atoms = [a.resi for a in gro.get_atoms(resn=probe_id)]
    resi_set = set(atoms)
    for resi in resi_set:
        mol = gro.get_atoms(resi=resi)
        com = center_of_mass(mol)
        virt_atom = GroAtom()
        virt_atom.point = com
        virt_atom.atomtype = "VIS"
        virt_atom.resn = mol[0].resn
        virt_atom.resi = resi
        gro.add_atom(virt_atom)

    logger.info(f"the system has {gro.molar(probe_id):.3f} M of {probe_id} cosolvents")
    return str(gro)


class VirtualAtomAdder:
    def __init__(self, sigma: float = 2, epsilon: float = 4.184e-6):
        # sigmaおよびepsilonは、VIS-VISのLJパラメータ。
        pass

    def add(self, system: System, target_probe: Probe) -> System:
        top_string = open(system.top.get()).read()
        gro_string = open(system.gro.get()).read()

        top_string_with_virt = addvirtatom2top(top_string, target_probe.cid)
        gro_string_with_virt = addvirtatom2gro(gro_string, target_probe.cid)

        return System.create_system_from_strings(top_string_with_virt, gro_string_with_virt)
