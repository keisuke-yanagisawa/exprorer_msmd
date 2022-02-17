#!/usr/bin/python3

import argparse
import configparser
import tempfile
from os.path import expanduser, expandvars, basename, dirname
import os
from subprocess import getoutput as gop

from scipy import constants
import jinja2

from utilities.util import expandpath


VERSION = "2.0.0"

TMP_PREFIX=".tmp"
EXT_FRCMOD=".frcmod"
EXT_PDB=".pdb"
EXT_INP=".in"

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

def protein_pdb_preparation(pdbfile):
    _, tmp1 = tempfile.mkstemp(prefix=TMP_PREFIX, suffix=EXT_PDB)
    _, tmp2 = tempfile.mkstemp(prefix=TMP_PREFIX, suffix=EXT_PDB)
    gop(f"grep -v OXT {pdbfile} > {tmp1}")
    gop(f"grep -v ANISOU {tmp1} > {tmp2}")
    return tmp2

def calculate_boxsize(pdbfile):
    tmp_prefix=f"{tmpdir}/{TMP_PREFIX}"
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


class Parmchk(object):
    def __init__(self, exe="parmchk2"):
        self.at_indices = {"gaff": 1, "gaff2": 2}
        self.exe = exe
    
    def set(self, mol2, at):
        self.at_id = self.at_indices[at]
        self.mol2 = mol2
        
    def run(self, frcmod=None):
        self.frcmod=frcmod
        if self.frcmod is None:
            _, self.frcmod = tempfile.mkstemp(prefix=TMP_PREFIX, suffix=EXT_FRCMOD)
        print(gop(f"{self.exe} -i {self.mol2} -f mol2 -o {frcmod} -s {self.at_id}"))

class Packmol(object):
    def __init__(self, exe="packmol"):
        self.exe = exe
        None
    def set(self, protein_pdb, cosolv_pdbs, box_size, molar):
        self.protein_pdb = protein_pdb
        self.cosolv_pdbs = cosolv_pdbs
        self.box_size = box_size
        self.molar = molar
    def run(self, box_pdb=None, seed=-1):
        self.box_pdb = box_pdb if not box_pdb is None \
                               else tempfile.mkstemp(prefix=TMP_PREFIX, suffix=EXT_PDB)[1]
        self.seed = seed

        # shorten path length to pdb file
        # too long path cannot be treated by packmol
        _, tmp_prot_pdb =  tempfile.mkstemp(prefix=TMP_PREFIX, suffix=EXT_PDB)
        print(gop(f"cp {self.protein_pdb} {tmp_prot_pdb}"))

        tmp_pdbs = [tempfile.mkstemp(prefix=TMP_PREFIX, suffix=EXT_PDB)[1] 
                    for _ in self.cosolv_pdbs]
        [print(gop(f"cp {src} {dst}")) for src, dst in zip(self.cosolv_pdbs, tmp_pdbs)]

        num = int(constants.N_A * self.molar * (self.box_size**3) * (10**-27))

        _, inp = tempfile.mkstemp(prefix=TMP_PREFIX, suffix=EXT_INP)
        with open(inp, "w") as fout:
            fout.write(TEMPLATE_PACKMOL_HEADER.format(output=self.box_pdb, prot=tmp_prot_pdb, seed=self.seed))
            fout.write("\n")
            for pdb in tmp_pdbs:
                fout.write(TEMPLATE_PACKMOL_STRUCT.format(cosolv=pdb, num=num, size=self.box_size/2))
                fout.write("\n")
        print(gop(f"{self.exe} < {inp}"))

    def __del__(self):
        os.remove(self.box_pdb)

class TLeap(object):
    def __init__(self, exe="tleap"):
        self.exe = exe
    def set(self, template_file, cids, cosolv_paths, frcmods, box_path, size, ssbonds, at):
        self.template_file = template_file
        self.cids = cids
        self.cosolv_paths = cosolv_paths
        self.frcmods = frcmods
        self.box_path = box_path
        self.size = size
        self.ssbonds = ssbonds
        self.at = at

    def run(self, oprefix):
        self.oprefix = oprefix
        self.parm7 = self.oprefix + ".parm7"
        self.rst7 = self.oprefix + ".rst7"

        _, inputfile = tempfile.mkstemp(prefix=TMP_PREFIX, suffix=EXT_INP)
        data = {
            "LIGAND_PARAM": f"leaprc.{self.at}",
            "SS_BONDS": zip(self.ssbonds[0::2], self.ssbonds[1::2]),
            "COSOLVENTS": [{"ID": cid, "PATH": cpath}
                        for cid, cpath in zip(self.cids, self.cosolv_paths)],
            "OUTPUT": self.oprefix,
            "SYSTEM_PATH": self.box_path,
            "COSOLV_FRCMODS": self.frcmods,
            "SIZE": self.size
        }

        env = jinja2.Environment(loader=jinja2.FileSystemLoader("/"))
        template = env.get_template(self.template_file)
        with open(inputfile, "w") as fout:
            fout.write(template.render(data))

        output = gop(f"{self.exe} -f {inputfile}")
        print(output)
        final_charge_info = [s.strip() for s in output.split("\n")
                            if s.strip().startswith("Total unperturbed charge")][0]
        self._final_charge_value = float(final_charge_info.split()[-1])

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
    parser.add_argument("--parmchk", dest="parmchk", default="parmchk2", help="path to parmchk")
    parser.add_argument("-tin", dest="tleap_input", required=True,
                        help="input file for tleap")
    parser.add_argument("-wat-ion-lst", dest="wat_ion_list", default="WAT,Na+,Cl-",
                        help="comma-separated water and ion list to be put on last in pdb entry")
    parser.add_argument("-no-rm-temp", action="store_true", dest="no_rm_temp_flag",
                        help="the flag not to remove all temporal files")
    parser.add_argument("-seed", default=-1, type=int)
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    # TODO: このファイルに対する相対パスでハードコーディングしてしまいたい。
    args.tleap_input = expandpath(args.tleap_input)


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
    cfrcmods = [f"{tmpdir}/.temp_cosolvent_{cid}.frcmod" for cid in cids]
    if "frcmod" in params["Cosolvent"] and params["Cosolvent"]["frcmod"] != "":
        cfrcmods = params["Cosolvent"]["frcmod"].split()
    packmol_obj = Packmol(args.packmol)
    packmol_obj.set(
        params["Protein"]["pdb"], 
        cpdbs,
        boxsize,
        float(params["Cosolvent"]["molar"])
    )
    while True:
        packmol_obj.run()

        if "frcmod" not in params["Cosolvent"] or params["Cosolvent"]["frcmod"] == "":
            for mol2, frcmod in zip(cmols, cfrcmods):
                parmchk = Parmchk(args.parmchk)
                parmchk.set(mol2, params["Cosolvent"]["atomtype"])
                parmchk.run(frcmod)

        # 2. amber tleap
        tleap_obj = TLeap(exe=args.tleap)
        tleap_obj.set(
            args.tleap_input, cids, cmols, cfrcmods,
            packmol_obj.box_pdb, boxsize, ssbonds,
            params["Cosolvent"]["atomtype"]
        )
        tleap_obj.run(args.output_prefix)
        system_charge = tleap_obj._final_charge_value

        if system_charge == 0:
            break
        else:
            print("the system is not neutral. generate system again")

    # 5. remove temporal files
    if not args.no_rm_temp_flag:
        print(gop(f"rm -r {tmpdir}"))
