import tempfile
import os
from subprocess import getoutput as gop
from scipy import constants
from utilities import const
from utilities.executable.execute import Command
import shutil
import jinja2

class Packmol(object):
    def __init__(self, exe="packmol", debug=False):
        self.exe = os.getenv("PACKMOL", "packmol")
        self.debug = debug

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

        shutil.copy2(self.protein_pdb, tmp_prot_pdb)

        tmp_pdbs = []
        for cosolv_pdb in self.cosolv_pdbs:
          _, tmp_pdb = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)
          shutil.copy2(cosolv_pdb, tmp_pdb)
          tmp_pdbs.append(tmp_pdb)

        num = int(constants.N_A * self.molar * (self.box_size**3) * (10**-27))

        _, inputfile = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)

        _, inp = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_INP)
        data = {
          "output":  self.box_pdb,
          "prot":    tmp_prot_pdb,
          "seed":    self.seed,
          "cosolvs": tmp_pdbs,
          "num":     num,
          "size":    self.box_size/2
        }

        env = jinja2.Environment(loader=jinja2.FileSystemLoader(f"{os.path.dirname(__file__)}/template"))
        template = env.get_template("packmol.in")
        with open(inp, "w") as fout:
            fout.write(template.render(data))
        command = Command(f"{self.exe} < {inp}")
        if self.debug:
          print(command)
        print(command.run())
        return self

    def __del__(self):
      if not self.debug:
        os.remove(self.box_pdb)