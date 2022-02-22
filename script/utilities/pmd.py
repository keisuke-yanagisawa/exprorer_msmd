import parmed as pmd
import tempfile

def _convert_top_only(intop, outtop):
  system = pmd.load_file(intop)
  system.save(outtop, overwrite=True)
  return outtop

def _convert_all(intop, outtop, inxyz, outxyz):
  system = pmd.load_file(intop, xyz=inxyz)
  system.save(outtop, overwrite=True)
  system.save(outxyz, overwrite=True)
  return outtop, outxyz

def convert(intop, outtop, inxyz=None, outxyz=None):
  if not inxyz is None:
    return _convert_all(intop, outtop, inxyz, outxyz)
  else:
    ret = _convert_top_only(intop, outtop)
    return ret, None