import tempfile
from typing import Optional, Tuple
from .variable import Path
import parmed as pmd


def _convert_top_only(intop: Path, outtop: Path) -> Path:
    system = pmd.load_file(intop.path)
    if outtop.path.endswith(".top"):  # この if ってどういう意味？
        protein_indices = [[list(t[1])[0] for t in system.split() if len(t[1]) == 1]]
        system.save(outtop.path, overwrite=True, combine=protein_indices)
    else:
        system.save(outtop.path, overwrite=True)
    return outtop


def _convert_all(intop: Path, outtop: Path, inxyz: Path, outxyz: Path):
    system = pmd.load_file(intop.path, xyz=inxyz.path)
    if outtop.path.endswith(".top"):  # この if ってどういう意味？
        # for dealing with multimer system
        # posre generator does not support multimer system
        # and thus we need to merge multimer into a record
        # Assumption: There is single molecule consisted with multiple residues -> proteins
        protein_indices = [[list(indices)[0] for parm, indices in system.split() if (len(indices) == 1 and len(parm.residues) > 1)]]
        system.save(outtop.path, overwrite=True, combine=protein_indices)
        system.save(outxyz.path, overwrite=True, combine=protein_indices)
    else:
        system.save(outtop.path, overwrite=True)
        system.save(outxyz.path, overwrite=True)
    return outtop, outxyz


def convert(intop: Path, inxyz: Optional[Path] = None) -> Tuple[Path, Optional[Path]]:
    outtop = Path(tempfile.mkstemp(suffix=".top")[1])
    outxyz = Path(tempfile.mkstemp(suffix=".gro")[1])
    if inxyz is not None:
        return _convert_all(intop, outtop, inxyz, outxyz)
    else:
        ret = _convert_top_only(intop, outtop)
        return ret, None
