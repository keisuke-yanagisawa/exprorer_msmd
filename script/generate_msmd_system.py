#!/usr/bin/python3

import tempfile
import os
from subprocess import getoutput as gop
from typing import Tuple

from script.utilities.executable import Parmchk, Packmol, TLeap
from script.utilities import const
from script.utilities.logger import logger

VERSION = "2.0.0"

tmp_leap = """
source leaprc.protein.ff14SB
source leaprc.water.tip3p

prot = loadPDB {pdbfile}

addIons2 prot Na+ 0
addIons2 prot Cl- 0
solvateBox prot TIP3PBOX 10.0

center prot
charge prot

saveAmberParm prot {tmp_prefix}.parm7 {tmp_prefix}.rst7
quit
"""


def protein_pdb_preparation(pdbfile: str) -> str:
    """
    remove OXT and ANISOU from pdbfile to avoid tleap error
    -----
    input
        pdbfile: path to pdb file
    output
        tmp1: path to pdb file without OXT and ANISOU
    """
    _, tmp1 = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)
    gop(f"grep -v OXT {pdbfile} | grep -v ANISOU > {tmp1}")
    return tmp1


def __calculate_boxsize(pdbfile: str) -> float:
    tmpdir = tempfile.mkdtemp()
    tmp_prefix = f"{tmpdir}/{const.TMP_PREFIX}"
    with open(f"{tmp_prefix}.in", "w") as fout:
        fout.write(tmp_leap.format(pdbfile=pdbfile, tmp_prefix=tmp_prefix))
    logger.info(gop(f"tleap -f {tmp_prefix}.in | tee {tmp_prefix}.in.result"))

    try:
        size = calculate_boxsize(f"{tmp_prefix}.rst7")
    except ValueError as e:
        logger.error("====tleap input commands====")
        logger.error(gop(f"cat {tmp_prefix}.in"))
        logger.error("====tleap output====")
        logger.error(gop(f"cat {tmp_prefix}.in.result"))
        raise e
    return size


def calculate_boxsize(rst7: str) -> float:
    """
    get longest box size from rst7 file
    """

    box_size_str = gop(f"tail -n 1 {rst7} | cut -c -36")
    box_size = [float(s) for s in box_size_str.split()]
    box_size = max(box_size)
    return box_size


def _create_frcmod(mol2file: str,
                   atomtype: str,
                   debug: bool = False) -> str:
    """
    create frcmod file from mol2 file
    -----
    input
        mol2file: path to mol2 file
        atomtype: atom type (GAFF / GAFF2)
    output
        cfrcmod: path to frcmod file
    """
    _, cfrcmod = tempfile.mkstemp(suffix=".frcmod")
    Parmchk(debug=debug) \
        .set(mol2file, atomtype) \
        .run(frcmod=cfrcmod)
    return cfrcmod


def create_system(setting_protein: dict,
                  setting_probe: dict,
                  probe_frcmod: str,
                  debug: bool = False,
                  seed: int = -1) -> Tuple[str, str]:
    """
    create system from protein and probe
    """
    pdbpath = protein_pdb_preparation(setting_protein["pdb"])
    boxsize = __calculate_boxsize(pdbpath)
    ssbonds = setting_protein["ssbond"]
    cmol = setting_probe["mol2"]
    cpdb = setting_probe["pdb"]
    cid = setting_probe["cid"]
    atomtype = setting_probe["atomtype"]
    probemolar = float(setting_probe["molar"])

    _, box_pdb = tempfile.mkstemp(suffix=".pdb")
    Packmol(debug=debug) \
        .set(pdbpath, cpdb, boxsize, probemolar) \
        .run(box_pdb, seed=seed)

    tleap_obj = TLeap(debug=debug) \
        .set(cid, cmol, probe_frcmod, box_pdb,
             boxsize, ssbonds, atomtype)

    while True:
        _, fileprefix = tempfile.mkstemp(suffix="")
        os.system(f"rm {fileprefix}")
        tleap_obj.run(fileprefix)
        system_charge = tleap_obj._final_charge_value

        if system_charge == 0:
            break
        else:
            logger.warn("the system is not neutral. generate system again")

    return tleap_obj.parm7, tleap_obj.rst7


def generate_msmd_system(setting: dict,
                         debug: bool = False,
                         seed: int = -1) -> Tuple[str, str]:
    """
    generate msmd system
    -----
    input
        setting: setting json (dict)
        debug: debug mode
        seed: random seed
    output
        parm7: path to parm7 file
        rst7: path to rst7 file
    """
    cfrcmod = _create_frcmod(setting["input"]["probe"]["mol2"],
                             setting["input"]["probe"]["atomtype"],
                             debug=debug)
    parm7, rst7 = create_system(setting["input"]["protein"],
                                setting["input"]["probe"],
                                cfrcmod, debug=debug,
                                seed=seed)
    return parm7, rst7
