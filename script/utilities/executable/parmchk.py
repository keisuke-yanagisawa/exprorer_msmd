import tempfile
import os
from utilities.executable.execute import Command
from utilities import const

class Parmchk(object):
  def __init__(self, debug=False):
    self.at_indices = {"gaff": 1, "gaff2": 2}
    self.exe = os.getenv("PARMCHK", "parmchk2")
    self.debug = debug
  
  def set(self, mol2, at):
    self.at_id = self.at_indices[at]
    self.mol2 = mol2
    return self
      
  def run(self, frcmod=None):
    self.frcmod=frcmod
    self.frcmod_user_define = self.frcmod is not None
    if self.frcmod is None:
      _, self.frcmod = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_FRCMOD)
    command = Command(f"{self.exe} -i {self.mol2} -f mol2 -o {self.frcmod} -s {self.at_id}")
    if self.debug:
      print(command)
    print(command.run())
    return self

  def __del__(self):
    if self.debug:
      print(self.frcmod)
    elif not self.frcmod_user_define:
      os.remove(self.frcmod)