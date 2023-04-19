#!/usr/bin/python3

import argparse
import tempfile
import os
from subprocess import getoutput as gop

from utilities import util
from utilities.executable import Parmchk, Packmol, TLeap
from utilities import const
from utilities.pmd import convert as pmd_convert
from utilities.logger import logger

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


def protein_pdb_preparation(pdbfile):
    _, tmp1 = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)
    gop(f"grep -v OXT {pdbfile} | grep -v ANISOU > {tmp1}")
    return tmp1


def __calculate_boxsize(pdbfile):
    tmpdir = tempfile.mkdtemp()
    tmp_prefix = f"{tmpdir}/{const.TMP_PREFIX}"
    with open(f"{tmp_prefix}.in", "w") as fout:
        fout.write(tmp_leap.format(pdbfile=pdbfile, tmp_prefix=tmp_prefix))
    logger.info(gop(f"tleap -f {tmp_prefix}.in | tee {tmp_prefix}.in.result"))
    return calculate_boxsize(f"{tmp_prefix}.rst7")


def calculate_boxsize(rst7):
    "get longest box size from rst7 file"
    box_size_str = gop(f"tail -n 1 {rst7} | cut -c -36")
    try:
        box_size = [float(s) for s in box_size_str.split()]
    except ValueError as e:
        logger.error(e)
        logger.error("cat leap.log")
        logger.error(open("leap.log").read())
        raise e
    box_size = max(box_size)
    return box_size


def create_frcmod(setting_probe, debug=False):
    cmol = setting_probe["mol2"]
    atomtype = setting_probe["atomtype"]
    _, cfrcmod = tempfile.mkstemp(suffix=".frcmod")
    Parmchk(debug=debug) \
        .set(cmol, atomtype) \
        .run(frcmod=cfrcmod)
    return cfrcmod


def create_system(setting_protein, setting_probe, probe_frcmod, debug=False, seed=-1):
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


def generate_msmd_system(setting, debug=False, seed=-1):
    cfrcmod = create_frcmod(setting["input"]["probe"], debug=debug)
    parm7, rst7 = create_system(setting["input"]["protein"],
                                setting["input"]["probe"],
                                cfrcmod, debug=debug,
                                seed=seed)
    return parm7, rst7


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="generate .parm7 and .rst7 files with cosolvent box")
    parser.add_argument("-setting-yaml", help="yaml file for the MSMD")
    parser.add_argument("-oprefix", dest="output_prefix", required=True)

    parser.add_argument("--seed", default=-1, type=int)
    parser.add_argument("--no-rm-temp", action="store_true", dest="no_rm_temp_flag",
                        help="the flag not to remove all temporal files")

    parser.add_argument("-v,--verbose", dest="verbose", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    if args.debug:
        logger.setLevel("debug")
    elif args.verbose:
        logger.setLevel("info")
    # else: logger level is "warn"

    setting = util.parse_yaml(args.setting_yaml)
    parm7, rst7 = generate_msmd_system(setting, args.debug, args.seed)
    pmd_convert(parm7, f"{args.output_prefix}.top",
                inxyz=rst7, outxyz=f"{args.output_prefix}.gro")
