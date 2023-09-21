import tempfile
from typing import Tuple
from .variable import Path
import parmed as pmd


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


def convert(intop: Path, inxyz: Path) -> Tuple[Path, Path]:
    outtop = Path(tempfile.mkstemp(suffix=".top")[1])
    outxyz = Path(tempfile.mkstemp(suffix=".gro")[1])
    return _convert_all(intop, outtop, inxyz, outxyz)
