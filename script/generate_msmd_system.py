#!/usr/bin/python3

import argparse
import configparser
import tempfile
from os.path import basename, dirname
import os
from subprocess import getoutput as gop

from scipy import constants
import jinja2
import parmed as pmd

from utilities.util import expandpath
from utilities.executable import Parmchk, Packmol, TLeap
from utilities import const


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
    _, tmp2 = tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)
    gop(f"grep -v OXT {pdbfile} > {tmp1}")
    gop(f"grep -v ANISOU {tmp1} > {tmp2}")
    return tmp2

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
    parser.add_argument("-prot_param", required=True,
                        help="parameter file of protein")
    parser.add_argument("-cosolv_param", required=True,
                        help="parameter file of cosolvent")
    parser.add_argument("-oprefix", dest="output_prefix", required=True)

    parser.add_argument("--seed", default=-1, type=int)
    parser.add_argument("--no-rm-temp", action="store_true", dest="no_rm_temp_flag",
                        help="the flag not to remove all temporal files")
    parser.add_argument("--wat-ion-lst", dest="wat_ion_list", default="WAT,Na+,Cl-,CA,MG,ZN,CU",
                        help="comma-separated water and ion list to be put on last in pdb entry")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    params = configparser.ConfigParser()
    params.read(expandpath(args.prot_param), "UTF-8")
    params.read(expandpath(args.cosolv_param), "UTF-8")
    #print(expandpath(args.prot_param), dir(params))
    params["Protein"]["pdb"]    = dirname(expandpath(args.prot_param))+"/"+basename(params["Protein"]["pdb"])
    params["Cosolvent"]["mol2"] = dirname(expandpath(args.cosolv_param))+"/"+basename(params["Cosolvent"]["mol2"])
    params["Cosolvent"]["pdb"]  = dirname(expandpath(args.cosolv_param))+"/"+basename(params["Cosolvent"]["pdb"])

    params["Protein"]["pdb"] = protein_pdb_preparation(params["Protein"]["pdb"])


    ssbonds = params["Protein"]["ssbond"].split()
    cmols = params["Cosolvent"]["mol2"].split()
    cpdbs = params["Cosolvent"]["pdb"].split()
    cids = params["Cosolvent"]["cid"].split()

    boxsize = calculate_boxsize(params["Protein"]["pdb"])


    # 0. preparation
    cosolv_box_size = ((constants.N_A * float(params["Cosolvent"]["molar"]))
                       * 1e-27)**(-1/3.0)  # /L -> /A^3 : 1e-27

    # 1. generate cosolvent box with packmol
    cfrcmods = []
    if "frcmod" in params["Cosolvent"] and params["Cosolvent"]["frcmod"] != "":
        cfrcmods = params["Cosolvent"]["frcmod"].split()
    else:
        for mol2 in cmols:
            parmchk = Parmchk(debug=args.debug)
            parmchk.set(mol2, params["Cosolvent"]["atomtype"]).run()
            cfrcmods.append(parmchk.frcmod)

    packmol_obj = Packmol(debug=args.debug)
    packmol_obj.set(
        params["Protein"]["pdb"], 
        cpdbs, boxsize,
        float(params["Cosolvent"]["molar"])
    )

    tleap_obj = TLeap(debug=args.debug)
    while True:
        packmol_obj.run()

        # 2. amber tleap
        tleap_obj.set(
            cids, cmols, cfrcmods,
            packmol_obj.box_pdb, boxsize, ssbonds,
            params["Cosolvent"]["atomtype"]
        ).run(args.output_prefix)
        system_charge = tleap_obj._final_charge_value

        if system_charge == 0:
            break
        else:
            print("the system is not neutral. generate system again")

    msmd_system = pmd.load_file(tleap_obj.parm7, xyz=tleap_obj.rst7)
    msmd_system.save(f"{args.output_prefix}.top", overwrite=True)
    msmd_system.save(f"{args.output_prefix}.gro", overwrite=True)

    # 5. remove temporal files
    if not args.no_rm_temp_flag:
        None # TODO: