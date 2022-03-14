#!/usr/bin/python3

import argparse
import configparser
import tempfile
from os.path import basename, dirname
import os
from subprocess import getoutput as gop

from scipy import constants
import jinja2

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

def calculate_boxsize(pdbfile):
    tmpdir = tempfile.mkdtemp()
    tmp_prefix=f"{tmpdir}/{const.TMP_PREFIX}"
    with open(f"{tmp_prefix}.in", "w") as fout:
        fout.write(tmp_leap.format(pdbfile=pdbfile, tmp_prefix=tmp_prefix))
    print(gop(f"tleap -f {tmp_prefix}.in | tee {tmp_prefix}.in.result"))
    center_str = gop(f"cat {tmp_prefix}.in.result | grep The\ center | cut -d: -f2 | sed s/,//g")
    center = [float(s) for s in center_str.split()]
    box_size_str = gop(f"tail -n 1 {tmp_prefix}.rst7 | cut -c -36")
    try:
        box_size = [float(s) for s in box_size_str.split()]
    except ValueError as e:
        print(e)
        print("cat leap.log")
        print(gop("cat leap.log"))
        exit(1)
        
    box_size = max(box_size)

    return box_size

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

    pdbpath = protein_pdb_preparation(setting["input"]["protein"]["pdb"])
    # probemolar = setting["input"]["probe"]["molar"]

    ssbonds    = setting["input"]["protein"]["ssbond"]
    cmol       = setting["input"]["probe"]["mol2"]
    cpdb       = setting["input"]["probe"]["pdb"]
    cid        = setting["input"]["probe"]["cid"]
    atomtype   = setting["input"]["probe"]["atomtype"]
    probemolar = float( setting["input"]["probe"]["molar"] )

    boxsize = calculate_boxsize(pdbpath)


    # 0. preparation
    cosolv_box_size = ((constants.N_A * float(probemolar))
                       * 1e-27)**(-1/3.0)  # /L -> /A^3 : 1e-27

    # 1. generate cosolvent box with packmol
    # input: 
    #   frcmod
    _, cfrcmod = tempfile.mkstemp(suffix=".frcmod")
    Parmchk(debug=args.debug) \
        .set(cmol, atomtype) \
        .run(frcmod=cfrcmod)
    
    _, box_pdb = tempfile.mkstemp(suffix=".pdb")
    Packmol(debug=args.debug) \
        .set(pdbpath, cpdb, boxsize, probemolar) \
        .run(box_pdb)

    tleap_obj = TLeap(debug=args.debug) \
                    .set(cid, cmol, cfrcmod, box_pdb, 
                         boxsize, ssbonds, atomtype)

    while True:
        tleap_obj.run(args.output_prefix)
        system_charge = tleap_obj._final_charge_value

        if system_charge == 0:
            break
        else:
            logger.warn("the system is not neutral. generate system again")

    pmd_convert(tleap_obj.parm7, f"{args.output_prefix}.top",
                inxyz=tleap_obj.rst7, outxyz=f"{args.output_prefix}.gro")

    # 5. remove temporal files
    if not args.no_rm_temp_flag:
        None # TODO: