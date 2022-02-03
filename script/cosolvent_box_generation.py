#!/usr/bin/python3

import argparse
import configparser
import tempfile
from os.path import expanduser, expandvars
from subprocess import getoutput as gop

from scipy import constants
import jinja2

from utilities.util import expandpath


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

TEMPLATE_PACKMOL_HEADER = """
seed {seed}
tolerance 2.0
output {output}
add_amber_ter
filetype pdb
structure {prot}
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  centerofmass
end structure
"""

TEMPLATE_PACKMOL_STRUCT = """
structure {cosolv}
  number {num}
  inside box -{size} -{size} -{size} {size} {size} {size}
end structure
"""

tmpdir = tempfile.mkdtemp()

def calculate_boxsize(pdbfile, tmp_prefix=".tmp"):
    tmp_prefix=f"{tmpdir}/{tmp_prefix}"
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


def run_parmchk(mol2, frcmod, at):
    at_id = {"gaff": 1, "gaff2": 2}[at]
    print(gop(f"parmchk2 -i {mol2} -f mol2 -o {frcmod} -s {at_id}"))


def gen_packmol_input(protein_pdb, cosolv_pdbs, box_pdb, inp, box_size, molar, seed=-1):
    # shorten path length to pdb file
    # too long path cannot be treated by packmol
    temp_protein_pdb = f"{tmpdir}/.temp_protein.pdb"
    print(gop(f"cp {protein_pdb} {temp_protein_pdb}"))

    temp_pdbs = [f"{tmpdir}/.temp_{i}.pdb" for i in range(len(cosolv_pdbs))]
    [gop(f"cp {src} {dst}") for src, dst in zip(cosolv_pdbs, temp_pdbs)]

    num = int(constants.N_A * molar * (box_size**3) * (10**-27))
    with open(inp, "w") as fout:
        fout.write(TEMPLATE_PACKMOL_HEADER.format(output=box_pdb, prot=temp_protein_pdb, seed=seed))
        fout.write("\n")
        for pdb in temp_pdbs:
            fout.write(TEMPLATE_PACKMOL_STRUCT.format(cosolv=pdb, num=num, size=box_size/2))
            fout.write("\n")


def run_packmol(packmol_path, inp):
    print(gop("%s < %s" % (packmol_path, inp)))


def gen_tleap_input(template_file, inputfile, cids, cosolv_paths, frcmods, box_path, size, oprefix, ssbonds, at):
    data = {
        "LIGAND_PARAM": f"leaprc.{at}",
        "SS_BONDS": zip(ssbonds[0::2], ssbonds[1::2]),
        "COSOLVENTS": [{"ID": cid, "PATH": cpath}
                       for cid, cpath in zip(cids, cosolv_paths)],
        "SIZE": size,
        "OUTPUT": oprefix,
        "SYSTEM_PATH": box_path,
        "COSOLV_FRCMODS": frcmods,
        "SIZE": size
    }

    env = jinja2.Environment(loader=jinja2.FileSystemLoader("/"))
    template = env.get_template(template_file)
    with open(inputfile, "w") as fout:
        fout.write(template.render(data))
    print(inputfile)

def run_tleap(tleap_path, inp):
    output = gop("%s -f %s" % (tleap_path, inp))
    print(output)
    final_charge_info = [s.strip() for s in output.split("\n")
                         if s.strip().startswith("Total unperturbed charge")][0]
    final_charge_value = float(final_charge_info.split()[-1])
    return final_charge_value

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="generate .parm7 and .rst7 files with cosolvent box")
    parser.add_argument("-prot_param", required=True,
                        help="parameter file of protein")
    parser.add_argument("-cosolv_param", required=True,
                        help="parameter file of cosolvent")
    parser.add_argument("-oprefix", dest="output_prefix", required=True)
    parser.add_argument("--packmol", dest="packmol", default="packmol", help="path to packmol")
    parser.add_argument("--tleap", dest="tleap", default="tleap", help="path to tleap")
    parser.add_argument("-tin", dest="tleap_input", required=True,
                        help="input file for tleap")
    parser.add_argument("-wat-ion-lst", dest="wat_ion_list", default="WAT,Na+,Cl-",
                        help="comma-separated water and ion list to be put on last in pdb entry")
    parser.add_argument("-no-rm-temp", action="store_true", dest="no_rm_temp_flag",
                        help="the flag not to remove all temporal files")
    parser.add_argument("-seed", default=-1, type=int)
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    params = configparser.ConfigParser()
    params.read(expandpath(args.prot_param), "UTF-8")
    params.read(expandpath(args.cosolv_param), "UTF-8")
    #print(expandpath(args.prot_param), dir(params))
    params["Protein"]["pdb"]    = expandpath(params["Protein"]["pdb"])
    params["Cosolvent"]["mol2"] = expandpath(params["Cosolvent"]["mol2"])
    params["Cosolvent"]["pdb"]  = expandpath(params["Cosolvent"]["pdb"])

    ssbonds = params["Protein"]["ssbond"].split()
    cmols = params["Cosolvent"]["mol2"].split()
    cpdbs = params["Cosolvent"]["pdb"].split()
    cids = params["Cosolvent"]["cid"].split()

    boxsize = calculate_boxsize(params["Protein"]["pdb"])


    # 0. preparation
    cosolv_box_size = ((constants.N_A * float(params["Cosolvent"]["molar"]))
                       * 1e-27)**(-1/3.0)  # /L -> /A^3 : 1e-27

    # 1. generate cosolvent box with packmol
    packmol_input = f"{tmpdir}/.temp_packmol.input"
    packmol_box_pdb = f"{tmpdir}/.temp_box.pdb"
    cfrcmods = [f"{tmpdir}/.temp_cosolvent_{cid}.frcmod" for cid in cids]
    if "frcmod" in params["Cosolvent"] and params["Cosolvent"]["frcmod"] != "":
        cfrcmods = params["Cosolvent"]["frcmod"].split()
    gen_packmol_input(params["Protein"]["pdb"], cpdbs,
                      packmol_box_pdb,
                      packmol_input,
                      boxsize,
                      float(params["Cosolvent"]["molar"]),
                      args.seed)
    while True:
        run_packmol(args.packmol, packmol_input)  # -> output: packmol_box_pdb
        temp_box = f"{tmpdir}/.temp_packmol_box_2.pdb"
        gop("grep -v OXT {} > {}".format(packmol_box_pdb, temp_box))

        if "frcmod" not in params["Cosolvent"] or params["Cosolvent"]["frcmod"] == "":
            for mol2, frcmod in zip(cmols, cfrcmods):
                run_parmchk(mol2, frcmod, params["Cosolvent"]["atomtype"])


        # 2. amber tleap
        tleap_input = f"{tmpdir}/.temp_tleap.in"
        gen_tleap_input(expandpath(args.tleap_input), tleap_input,
                        cids, cmols,
                        cfrcmods,
                        temp_box,
                        boxsize,
                        args.output_prefix,
                        ssbonds,
                        params["Cosolvent"]["atomtype"])
        system_charge = run_tleap(args.tleap, tleap_input)
        if system_charge == 0:
            break
        else:
            print("the system is not neutral. generate system again")

    # 5. remove temporal files
    if not args.no_rm_temp_flag:
        print(gop(f"rm -r {tmpdir}"))
