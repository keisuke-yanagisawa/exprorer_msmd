import tempfile
import os
from subprocess import getoutput as gop
from scipy import constants
from utilities import const

TEMPLATE_PACKMOL_HEADER = """
seed {seed}
tolerance 2.0
output {output}
add_amber_ter
filetype pdb
structure {prot}
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  centerofmass
end structure
"""

TEMPLATE_PACKMOL_STRUCT = """
structure {cosolv}
  number {num}
  inside box -{size} -{size} -{size} {size} {size} {size}
end structure
"""


class Packmol(object):
    def __init__(self, exe="packmol"):
        self.exe = exe
        None

    def set(self, protein_pdb, cosolv_pdbs, box_size, molar):
        self.protein_pdb = protein_pdb
        self.cosolv_pdbs = cosolv_pdbs
        self.box_size = box_size
        self.molar = molar
        return self

    def run(self, box_pdb=None, seed=-1):
        self.box_pdb = box_pdb if not box_pdb is None \
                               else tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)[1]
        self.seed = seed

        # shorten path length to pdb file
        # too long path cannot be treated by packmol
        _, tmp_prot_pdb =  tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)
        print(gop(f"cp {self.protein_pdb} {tmp_prot_pdb}"))

        tmp_pdbs = [tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)[1] 
                    for _ in self.cosolv_pdbs]
        [print(gop(f"cp {src} {dst}")) for src, dst in zip(self.cosolv_pdbs, tmp_pdbs)]

        num = int(constants.N_A * self.molar * (self.box_size**3) * (10**-27))

        _, inp = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)
        with open(inp, "w") as fout:
            fout.write(TEMPLATE_PACKMOL_HEADER.format(output=self.box_pdb, prot=tmp_prot_pdb, seed=self.seed))
            fout.write("\n")
            for pdb in tmp_pdbs:
                fout.write(TEMPLATE_PACKMOL_STRUCT.format(cosolv=pdb, num=num, size=self.box_size/2))
                fout.write("\n")
        print(gop(f"{self.exe} < {inp}"))
        return self

    def __del__(self):
        os.remove(self.box_pdb)