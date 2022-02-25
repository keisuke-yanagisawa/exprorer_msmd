#!/usr/bin/python3

import argparse
import configparser
import tempfile
from os.path import basename, dirname
import os
from subprocess import getoutput as gop

from scipy import constants
import jinja2

from utilities.pmd import convert as pmd_convert
from utilities.util import expandpath
from utilities.executable import Cpptraj
from utilities import const, GridUtil
from gridData import Grid
import numpy as np

VERSION = "1.0.0"

def check(gs):
    return True
def grid_max(gs):
    # TODO: check all grids have the same voxel fields
    if not check(gs):
        print("ERROR: Grid(s) have different voxel fields")
        exit(1)
    
    ret = gs[0]
    ret.grid = np.max([g.grid for g in gs], axis=0)
    return ret

def gen_max_pmap(inpaths, outpath):
    gs = [Grid(n) for n in inpaths]
    max_pmap = grid_max(gs)
    max_pmap.export(outpath, type="double")
    return outpath


DEBUG = None # global variable 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="generate PMAPs")
    parser.add_argument("pmap_paths", nargs="+",
                        help="PMAPs will be integrated")
    parser.add_argument("mappmap_path", 
                        help="destination path to output max-PMAP")

    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()
    DEBUG = args.debug # assign to a global variable 

    max_pmap_path = gen_max_pmap(args.pmap_paths, args.mappmap_path)

    # PMAPを各runに対して作成 - gen_pmap.sh
    #  この時点で当然probabilityになっているべき。
    # max-PMAPを作成 - max_pmap.sh
    # cpptrajでほげほげする - cpptraj.py
    # environmentを全列挙する - extract_environment.sh
    #   このあたりは計算時間がトンデモなくかかるので、1つのスクリプトにまとめない方が良い気がする
    # 構造の重ね合わせ - superimpose.sh
    # meshの描画 - make_meshes.py / environments_visualization.sh


    # parser.add_argument("-prot_param", required=True,
    #                     help="parameter file of protein")
    # parser.add_argument("-cosolv_param", required=True,
    #                     help="parameter file of cosolvent")
    # parser.add_argument("-oprefix", dest="output_prefix", required=True)

    # parser.add_argument("--seed", default=-1, type=int)
    # parser.add_argument("--no-rm-temp", action="store_true", dest="no_rm_temp_flag",
    #                     help="the flag not to remove all temporal files")
    # parser.add_argument("--wat-ion-lst", dest="wat_ion_list", default="WAT,Na+,Cl-,CA,MG,ZN,CU",
    #                     help="comma-separated water and ion list to be put on last in pdb entry")

# import argparse
# import configparser
# from subprocess import getoutput as gop

# import parmed as pmd

# from utilities import util


# def gen_parm7(in_top, in_gro, out_parm7):
#     gromacs = pmd.load_file(in_top, xyz=in_gro)
#     gromacs.save(out_parm7, overwrite=True)


# def do_cpptraj(top, traj, ref, cid, prefix, n=100, d=1, traj_start=1, traj_stop="last", traj_offset=1):
#     traj = util.expandpath(traj)
#     ref = util.expandpath(ref)
    
#     with open(".cpptraj.in", "w") as fout:
#         fout.write(f"""
# trajin {traj} {traj_start} {traj_stop} {traj_offset}
# parm {ref}
# reference {ref} parm {ref} [REF]

# unwrap @CA&(!:CA)&(!:{cid})
# center @CA&(!:CA)&(!:{cid})
# fiximagedbonds
# autoimage

# rms ToREF ref [REF] @CA&(!:CA)&(!:{cid}) @CA&(!:CA)&(!:{cid}) out rmsd.dat


# grid map_{prefix}_nV.dx {n} {d} {n} {d} {n} {d} :{cid}&(!@VIS)
# grid map_{prefix}_nVH.dx {n} {d} {n} {d} {n} {d} :{cid}&(!@VIS)&(!@H*)
# grid map_{prefix}_V.dx {n} {d} {n} {d} {n} {d} :{cid}@VIS
# grid map_{prefix}_O.dx {n} {d} {n} {d} {n} {d} :{cid}@O*

# trajout {prefix}_position_check.pdb start 1 stop 1 offset 1
# trajout {prefix}_position_check2.pdb offset 100
# trajout {prefix}_last_frame.pdb start 2001
# go
# exit
# """)
#     log = gop(f"cpptraj -p {top} < .cpptraj.in")
#     print(log)


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="gen probability map")
#     parser.add_argument("-p, --topology", dest="topology", required=True)
#     parser.add_argument("-y, --trajin", dest="trajin", required=True)
#     parser.add_argument("-c, --reference", dest="reference", required=True)
#     parser.add_argument("-o, --outprefix", dest="outprefix", default="map_1-10000")
#     parser.add_argument("--start", dest="start", default=1)
#     parser.add_argument("--stop", dest="stop", default="last")
#     parser.add_argument("--offset", dest="offset", default=1)
#     parser.add_argument("prot_param")
#     parser.add_argument("cosolv_param")

#     args = parser.parse_args()
#     params = configparser.ConfigParser()
#     params.read(args.prot_param)
#     params.read(args.cosolv_param)

#     temp_parm7 = ".temp.parm7"
#     gen_parm7(args.topology,
#               args.reference,
#               temp_parm7)

#     if "ReferenceStructure" not in params:
#         params["ReferenceStructure"] = params["Protein"]
#     do_cpptraj(temp_parm7,
#                args.trajin,
#                params["ReferenceStructure"]["pdb"],
#                params["Cosolvent"]["cid"],
#                args.outprefix,
#                traj_start=args.start,
#                traj_stop=args.stop,
#                traj_offset=args.offset)
