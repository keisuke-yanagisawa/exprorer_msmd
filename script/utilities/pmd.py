from pathlib import Path
from typing import Optional

import parmed as pmd


def _convert_top_only(intop: Path, outtop: Path) -> Path:
    system = pmd.load_file(str(intop))
    if outtop.suffix == ".top":
        protein_indices = [[list(t[1])[0] for t in system.split() if len(t[1]) == 1]]
        system.save(str(outtop), overwrite=True, combine=protein_indices)
    else:
        system.save(str(outtop), overwrite=True)
    return outtop


def _convert_all(intop: Path, outtop: Path, inxyz: Path, outxyz: Path) -> tuple[Path, Path]:
    system = pmd.load_file(str(intop), xyz=str(inxyz))
    if outtop.suffix == ".top":
        # for dealing with multimer system
        # posre generator does not support multimer system
        # and thus we need to merge multimer into a record
        # Assumption: There is single molecule consisted with multiple residues -> proteins
        protein_indices = [
            [list(indices)[0] for parm, indices in system.split() if (len(indices) == 1 and len(parm.residues) > 1)]
        ]
        system.save(str(outtop), overwrite=True, combine=protein_indices)
        system.save(str(outxyz), overwrite=True, combine=protein_indices)
    else:
        system.save(str(outtop), overwrite=True)
        system.save(str(outxyz), overwrite=True)
    return outtop, outxyz


def convert(intop: Path, outtop: Path, inxyz: Optional[Path] = None, outxyz: Optional[Path] = None):
    if inxyz is not None and outxyz is not None:
        return _convert_all(intop, outtop, inxyz, outxyz)
    else:
        ret = _convert_top_only(intop, outtop)
        return ret, None
