import parmed as pmd


def _convert_top_only(intop, outtop):
    system = pmd.load_file(intop)
    if outtop.endswith(".top"):
        protein_indices = [[list(t[1])[0] for t in system.split() if len(t[1]) == 1]]
        system.save(outtop, overwrite=True, combine=protein_indices)
    else:
        system.save(outtop, overwrite=True)
    return outtop


def _convert_all(intop, outtop, inxyz, outxyz):
    system = pmd.load_file(intop, xyz=inxyz)
    if outtop.endswith(".top"):
        # for dealing with multimer system
        # posre generator does not support multimer system
        # and thus we need to merge multimer into a record
        # Assumption: There is single molecule -> proteins
        #             There are multiple the same molecule -> probes and waters
        protein_indices = [[list(t[1])[0] for t in system.split() if len(t[1]) == 1]]
        system.save(outtop, overwrite=True, combine=protein_indices)
        system.save(outxyz, overwrite=True, combine=protein_indices)
    else:
        system.save(outtop, overwrite=True)
        system.save(outxyz, overwrite=True)
    return outtop, outxyz


def convert(intop, outtop, inxyz=None, outxyz=None):
    if inxyz is not None:
        return _convert_all(intop, outtop, inxyz, outxyz)
    else:
        ret = _convert_top_only(intop, outtop)
        return ret, None
