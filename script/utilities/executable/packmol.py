import collections
import tempfile
import os
from subprocess import getoutput as gop
from scipy import constants
from .. import const
from .execute import Command
import shutil
import jinja2
from ..logger import logger
from ..Bio.PDB import get_structure
from ..scipy.spatial import estimate_volume
import numpy as np

# obtained from biopython/Bio/PDB/SASA.py
__ATOMIC_RADII = collections.defaultdict(lambda: 2.0)
__ATOMIC_RADII.update(
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

    def set(self, protein_pdb: str, cosolv_pdb: str, box_size, molar: float):
        """
        protein_pdb: protein pdb file
        cosolv_pdb: cosolv pdb file
        """
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
      radii  = []
      for atom in get_structure(self.protein_pdb).get_atoms():
        if not atom.get_parent().resname in ["HOH", "WAT"]:
          coords.append(atom.get_coord())
          radii.append(__ATOMIC_RADII[atom.element])
      coords = np.array(coords)
      radii  = np.array(radii)
      radii += __ATOMIC_RADII["C"] # solvents' VdW radius: estimated by carbon radius.
      return estimate_volume(coords, radii)


    def run(self, box_pdb=None, seed=-1):
        self.box_pdb = box_pdb if not box_pdb is None \
                               else tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)[1]
        self.box_pdb_user_define = box_pdb is not None
        self.seed = seed

        # shorten path length to pdb file
        # too long path cannot be treated by packmol
        _, tmp_prot_pdb =  tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)

        shutil.copy2(self.protein_pdb, tmp_prot_pdb)

        _, tmp_pdb = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)
        shutil.copy2(self.cosolv_pdb, tmp_pdb)

        # determine the number of probe molecules
        # note: the generated system will shrink because of NPT-ensemble,
        #       so the number of probe molecules is decreased by the factor of 0.8.
        factor = 0.8
        protein_volume = self.__estimate_protein_exclute_volume()
        num = int(constants.N_A * self.molar * (self.box_size**3 - protein_volume) * (10**-27) * factor)

        _, inputfile = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)

        _, inp = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)
        data = {
          "output": self.box_pdb,
          "prot":   tmp_prot_pdb,
          "seed":   self.seed,
          "probe":  tmp_pdb,
          "num":    num,
          "size":   self.box_size/2
        }

        env = jinja2.Environment(loader=jinja2.FileSystemLoader(f"{os.path.dirname(__file__)}/template"))
        template = env.get_template("packmol.in")
        with open(inp, "w") as fout:
            fout.write(template.render(data))
        command = Command(f"{self.exe} < {inp}")
        logger.debug(command)
        logger.info(command.run())
        return self

    def __del__(self):
      if not self.debug:
        if not self.box_pdb_user_define:
          os.remove(self.box_pdb)