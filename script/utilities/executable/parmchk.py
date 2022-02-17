import tempfile
import os
from subprocess import getoutput as gop
from utilities import const

class Parmchk(object):
    def __init__(self, exe="parmchk2"):
        self.at_indices = {"gaff": 1, "gaff2": 2}
        self.exe = exe
    
    def set(self, mol2, at):
        self.at_id = self.at_indices[at]
        self.mol2 = mol2
        return self
        
    def run(self, frcmod=None):
        self.frcmod=frcmod
        if self.frcmod is None:
            _, self.frcmod = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_FRCMOD)
        print(gop(f"{self.exe} -i {self.mol2} -f mol2 -o {frcmod} -s {self.at_id}"))
        return self

    def __del__(self):
      os.remove(self.frcmod)