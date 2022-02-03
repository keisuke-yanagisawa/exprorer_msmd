import argparse
import configparser
from subprocess import getoutput as gop

import parmed as pmd

from utilities import util


def gen_parm7(in_top, in_gro, out_parm7):
    gromacs = pmd.load_file(in_top, xyz=in_gro)
    gromacs.save(out_parm7, overwrite=True)


def do_cpptraj(top, traj, ref, cid, prefix, n=100, d=1, traj_start=1, traj_stop="last", traj_offset=1):
    traj = util.expandpath(traj)
    ref = util.expandpath(ref)
    
    with open(".cpptraj.in", "w") as fout:
        fout.write(f"""
trajin {traj} {traj_start} {traj_stop} {traj_offset}
parm {ref}
reference {ref} parm {ref} [REF]

unwrap @CA&(!:CA)&(!:{cid})
center @CA&(!:CA)&(!:{cid})
fiximagedbonds
autoimage

rms ToREF ref [REF] @CA&(!:CA)&(!:{cid}) @CA&(!:CA)&(!:{cid}) out rmsd.dat


grid map_{prefix}_nV.dx {n} {d} {n} {d} {n} {d} :{cid}&(!@VIS)
grid map_{prefix}_nVH.dx {n} {d} {n} {d} {n} {d} :{cid}&(!@VIS)&(!@H*)
grid map_{prefix}_V.dx {n} {d} {n} {d} {n} {d} :{cid}@VIS
grid map_{prefix}_O.dx {n} {d} {n} {d} {n} {d} :{cid}@O*

trajout {prefix}_position_check.pdb start 1 stop 1 offset 1
trajout {prefix}_position_check2.pdb offset 100
trajout {prefix}_last_frame.pdb start 2001
go
exit
""")
    log = gop(f"cpptraj -p {top} < .cpptraj.in")
    print(log)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="gen probability map")
    parser.add_argument("-p, --topology", dest="topology", required=True)
    parser.add_argument("-y, --trajin", dest="trajin", required=True)
    parser.add_argument("-c, --reference", dest="reference", required=True)
    parser.add_argument("-o, --outprefix", dest="outprefix", default="map_1-10000")
    parser.add_argument("--start", dest="start", default=1)
    parser.add_argument("--stop", dest="stop", default="last")
    parser.add_argument("--offset", dest="offset", default=1)
    parser.add_argument("prot_param")
    parser.add_argument("cosolv_param")

    args = parser.parse_args()
    params = configparser.ConfigParser()
    params.read(args.prot_param)
    params.read(args.cosolv_param)

    temp_parm7 = ".temp.parm7"
    gen_parm7(args.topology,
              args.reference,
              temp_parm7)

    if "ReferenceStructure" not in params:
        params["ReferenceStructure"] = params["Protein"]
    do_cpptraj(temp_parm7,
               args.trajin,
               params["ReferenceStructure"]["pdb"],
               params["Cosolvent"]["cid"],
               args.outprefix,
               traj_start=args.start,
               traj_stop=args.stop,
               traj_offset=args.offset)
