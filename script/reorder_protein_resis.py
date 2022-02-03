import argparse
import configparser
import tempfile
import os

from Bio import PDB

from utilities.util import expandpath


def gen_resis_table(not_reordered, reordered):
    old_resis = []
    for res in not_reordered.get_residues():
        old_resis.append(res.id[1])

    new_resis = []
    for res in not_reordered.get_residues():
        new_resis.append(res.id[1])

    assert(len(old_resis) == len(new_resis))

    old_to_new = {}
    for old, new in zip(old_resis, new_resis):
        old_to_new[old] = new

    return old_to_new

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="generate residue number-reordered protein PDB and protein config file")
    parser.add_argument("-iconf", required=True)
    parser.add_argument("-oconf", required=True)
    parser.add_argument("-opdb",  required=True)
    parser.add_argument("--gromacs", default="gmx")
    parser.add_argument("--reduce", default="reduce",
                        help="reduce executable included in AmberTools")
    args = parser.parse_args()
    GMX=args.gromacs
    REDUCE=args.reduce

    params = configparser.ConfigParser()
    params.read(expandpath(args.iconf), "UTF-8")
    tmp = tempfile.mkstemp()[1]
    
    params["Protein"]["pdb"] = expandpath(params["Protein"]["pdb"])
    output_pdb_filepath = args.opdb
    os.system(f"{REDUCE} -Trim {params['Protein']['pdb']} > {tmp}.noH.pdb")
    os.system(f"""{GMX} trjconv -f {tmp}.noH.pdb -s {tmp}.noH.pdb -o {tmp}.center.pdb -center -boxcenter zero << EOF
Protein
System
EOF""")
    os.system(f"{GMX} editconf -f {tmp}.center.pdb -o {output_pdb_filepath} -resnr 1")

    old_pdb = PDB.PDBParser().get_structure("old", params["Protein"]["pdb"])
    new_pdb = PDB.PDBParser().get_structure("new", output_pdb_filepath)
    old_to_new = gen_resis_table(old_pdb, new_pdb)

    length = len([res for res in new_pdb.get_residues() if res.id[0] == " "])
    params["Protein"]["resi_st"] = "1"
    params["Protein"]["resi_ed"] = f"{length}"
    params["Protein"]["pdb"] = output_pdb_filepath
    new_ssbonds = [old_to_new[int(s)] for s in params["Protein"]["ssbond"].split()]
    params["Protein"]["ssbond"] = " ".join([str(s) for s in new_ssbonds])
    new_binding_residues = [old_to_new[int(s)] for s in params["Protein"]["binding_site_residues"].split()]
    params["Protein"]["binding_site_residues"] = " ".join([str(s) for s in new_binding_residues])

    with open(args.oconf, "w") as fout:
        params.write(fout)
