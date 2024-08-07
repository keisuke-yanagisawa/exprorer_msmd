import collections
import os
import shutil
import tempfile
import warnings
from pathlib import Path
from typing import Optional

import jinja2
import numpy as np
from scipy import constants

from .. import const
from ..Bio.PDB import get_structure
from ..logger import logger
from ..scipy.spatial_func import estimate_volume
from .execute import Command

# obtained from biopython/Bio/PDB/SASA.py
_ATOMIC_RADII = collections.defaultdict(lambda: 2.0)
_ATOMIC_RADII.update(
    {
        "H": 1.200,
        "HE": 1.400,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
    }
)


class Packmol(object):
    def __init__(self, exe="packmol", debug=False):
        self.exe = os.getenv("PACKMOL", "packmol")
        self.debug = debug
        self.box_pdb = None

    def __is_cosolvent_pdb(self, pdb: Path) -> bool:
        """
        check if pdb is cosolvent pdb
        Assumption: cosolvent pdb has only one residue
        """
        struct = get_structure(pdb)
        residues = [residue for residue in struct.get_residues()]
        return len(residues) == 1

    def __is_protein_pdb(self, pdb: Path) -> bool:
        """
        check if pdb is protein pdb
        Assumption: protein has more than one residue
        """
        struct = get_structure(pdb)
        residues = [residue for residue in struct.get_residues()]
        return len(residues) > 1

    def set(self, protein_pdb: Path, cosolv_pdb: Path, box_size: float, molar: float):
        """
        protein_pdb: protein pdb file
        cosolv_pdb: cosolv pdb file
        """
        if molar < 0:
            raise ValueError("molar must be zero or positive value")
        elif molar == 0:
            warnings.warn("molar is 0.0. No cosolvent will be added.", RuntimeWarning)
        elif molar > 55:
            warnings.warn(
                "given molar is greater than 55.0 M but truncated to 55.0 M. The concentration is extremely higher than the usual setting in MSMD simulation.",
                RuntimeWarning,
            )
            molar = 55.0
        elif molar > 1.0:
            warnings.warn(
                "molar is greater than 1.0 M. The concentration is higher than the usual setting in MSMD simulation.",
                RuntimeWarning,
            )
        if os.path.splitext(protein_pdb)[1] != ".pdb":
            raise ValueError(f"protein_pdb {protein_pdb} must be pdb file")
        if os.path.splitext(cosolv_pdb)[1] != ".pdb":
            raise ValueError(f"cosolv_pdb {cosolv_pdb} must be pdb file")
        if not self.__is_cosolvent_pdb(cosolv_pdb):
            raise RuntimeError(f"cosolv_pdb {cosolv_pdb} must have only one residue")
        if not self.__is_protein_pdb(protein_pdb):
            raise RuntimeError(f"protein_pdb {protein_pdb} must have more than one residue")
        self.protein_pdb = protein_pdb
        self.cosolv_pdb = cosolv_pdb
        self.box_size = box_size
        self.molar = molar
        return self

    def __estimate_protein_exclute_volume(self) -> float:
        """
        VdW半径に基づいてタンパク質の排除体積を計算
        タンパク質は原子種類に応じて処理し、
        溶媒は炭素原子1つ分の大きさであるとする
        """
        coords = []
        radii = []
        for atom in get_structure(self.protein_pdb).get_atoms():
            if not atom.get_parent().resname in ["HOH", "WAT"]:
                coords.append(atom.get_coord())
                radii.append(_ATOMIC_RADII[atom.element])
        coords = np.array(coords)
        radii = np.array(radii)
        radii += _ATOMIC_RADII["C"]  # solvents' VdW radius: estimated by carbon radius.
        return estimate_volume(coords, radii)

    def run(self, box_pdb: Optional[Path] = None, seed=-1):
        self.box_pdb = (
            box_pdb if box_pdb is not None else Path(tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)[1])
        )
        self.seed = seed

        # shorten path length to pdb file
        # too long path cannot be treated by packmol
        tmp_prot_pdb = Path(tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)[1])

        shutil.copy2(self.protein_pdb, tmp_prot_pdb)

        tmp_pdb = Path(tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)[1])
        shutil.copy2(self.cosolv_pdb, tmp_pdb)

        # determine the number of probe molecules
        # note: the generated system will shrink because of NPT-ensemble,
        #       so the number of probe molecules is decreased by the factor of 0.8.
        factor = 0.8
        protein_volume = self.__estimate_protein_exclute_volume()
        num = int(constants.N_A * self.molar * (self.box_size**3 - protein_volume) * (10**-27) * factor)

        _, inputfile = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)

        _, inp = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)
        probes = [{"pdb": tmp_pdb, "num": num, "size": self.box_size / 2}]

        data = {
            "output": self.box_pdb,
            "prot": tmp_prot_pdb,
            "seed": self.seed,
            "probes": [probe for probe in probes if probe["num"] > 0],
        }

        env = jinja2.Environment(loader=jinja2.FileSystemLoader(f"{os.path.dirname(__file__)}/template"))
        template = env.get_template("packmol.in")
        with open(inp, "w") as fout:
            fout.write(template.render(data))
        if self.debug:
            logfile = "packmol.log"
        else:
            _, logfile = tempfile.mkstemp(suffix=".log")
        command = Command(f"{self.exe} < {inp} > {logfile}")
        logger.debug(command)
        command.run()

        stdout_str = open(logfile).read()
        if stdout_str.find("ENDED WITHOUT PERFECT PACKING") != -1:
            raise RuntimeError(f"packmol failed: \n{stdout_str}")
        elif stdout_str.find("There are only fixed molecules") != -1:
            warnings.warn("There is only a protein molecule.", RuntimeWarning)
        logger.info(stdout_str)
        return self
