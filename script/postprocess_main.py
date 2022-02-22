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
import gridData
import numpy as np

VERSION = "1.0.0"

def mask_generator(ref_struct, reference_grid, distance=None):
    mask = GridUtil.gen_distance_grid(reference_grid, ref_struct)
    # print(np.max(mask.grid), np.min(mask.grid), distance)
    if distance is not None:
        mask.grid = mask.grid < distance
    else:
        mask.grid = mask.grid < np.inf
    return mask

def convert_to_proba(g, mask_grid=None):
    if mask_grid is not None:
        values = g.grid[np.where(mask_grid)]
        # print(np.sum(g.grid), np.sum(values), np.where(mask_grid))
        values /= np.sum(values)
        g.grid = np.full_like(g.grid, np.min([np.min(values), -1])) #TODO: I could not remember why it is needed...??
        g.grid[np.where(mask_grid)] = values
    else:
        g.grid /= np.sum(g.grid)
    return g

def convert_to_pmap(grid_path, ref_struct, valid_distance):
    grid = gridData.Grid(grid_path)
    mask = mask_generator(ref_struct, grid, valid_distance)
    pmap = convert_to_proba(grid, mask.grid)
    
    pmap_path = grid_path + "_pmap.dx" # TODO: terrible naming
    pmap.export(pmap_path, type="double")
    return pmap_path


def run_gen_pmap(basedir, prot_param, cosolv_param, valid_distance, debug=False):
    params = configparser.ConfigParser()
    params.read(expandpath(prot_param), "UTF-8")
    params.read(expandpath(cosolv_param), "UTF-8")
    if "ReferenceStructure" not in params:
        params["ReferenceStructure"] = params["Protein"]

    trajectory = expandpath(basedir) + "/" + "simulation" + "/" + "TEST_PROJECT.xtc" # TODO: TEST_PROJECT should be "PREFIX"
    topology   = expandpath(basedir) + "/" + "top" + "/" + "TEST_PROJECT.top" # TODO: TEST_PROJECT should be "PREFIX"
    ref_struct = dirname(expandpath(prot_param))+"/"+basename(params["ReferenceStructure"]["pdb"])
    probe_id   = params["Cosolvent"]["cid"]

    cpptraj_obj = Cpptraj(debug=debug)
    cpptraj_obj.set(topology, trajectory, ref_struct, probe_id)
    cpptraj_obj.run(
        basedir=basedir, 
        prefix="test"
    )
    
    pmap_paths = []
    for grid_path in cpptraj_obj.grids:
        pmap_path = convert_to_pmap(grid_path, ref_struct, valid_distance)
        pmap_paths.append(pmap_path)
    
    return pmap_paths



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="postprocessing")
    parser.add_argument("-basedir", required=True,
                        help="objective directory")
    parser.add_argument("prot_param")
    parser.add_argument("cosolv_param")

    parser.add_argument("-d,--distance-threshold", dest="d", metavar="d", default=5, type=int,
                        help="distance from protein atoms.")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args()

    pmap_paths = run_gen_pmap(args.basedir, args.prot_param, args.cosolv_param, args.d, args.debug)

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
