#! /usr/bin/python3

import argparse
import configparser
import pathlib
import subprocess
import shutil
import jinja2

VERSION = "2.1.0"

template_minimize = """
define       = {defines}
;
cutoff-scheme = Verlet
integrator   = steep
nsteps       = {steps}
emstep       = {emstep}
;
nstlog       = 1
;
nstlist      = 100
ns_type      = grid
pbc          = xyz
rlist        = 1.0
;
coulombtype  = PME
rcoulomb     = 1.0
;
vdwtype      = cut-off
rvdw         = 1.0
;
optimize_fft = yes
;
constraints  = none
"""

template_heat = """
; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)
define       = %s
; RUN CONTROL PARAMETERS
integrator   = md
dt           = %s
nsteps       = %s
;
nstxtcout     = 500
xtc-precision = 10000
nstlog        = 500
nstenergy     = 500
;
nstlist      = 5
ns_type      = grid
pbc          = xyz
rlist        = 0.8
;
coulombtype  = PME
rcoulomb     = 0.8
;
vdwtype      = cut-off
rvdw         = 0.8
dispcorr     = enerpres
;
optimize_fft = yes
;
tcoupl       = v-rescale
tc_grps      = Non-Water Water
tau_t        = 0.1 0.1
ref_t        = 300.0 300.0
;
pcoupl       = %s
pcoupltype   = isotropic
tau_p        = 1.0
compressibility = 4.5e-5
ref_p        = 1.0
refcoord_scaling = all
;
%s
;
continuation = %s ; Restarting after previous RUN
constraints  = hbonds
constraint_algorithm = LINCS
%s
"""

template_pr = """
define            = {defines}
; Run parameters
integrator      = md            ; leap-frog integrator
dt              = 0.002     ; 2 fs
nsteps          = {nsteps}
cutoff-scheme   = Verlet

; Output control
nstxtcout       = {snapshot_interval}          ; save coordinates every N ps
xtc-precision   = 10000
nstenergy       = 2500              ; save energies every 5 ps
nstlog          = 5000        ; update log file every 10 ps

; Bond parameters
continuation    = yes           ; Restarting after previous RUN
constraint_algorithm = LINCS    ; holonomic constraints
constraints          = hbonds  ; Convert the bonds with H-atoms to constraints
lincs_iter           = 1            ; accuracy of LINCS
lincs_order          = 4              ; also related to accuracy

; Langevin Dynamics
bd_fric = 2.0           ; Brownian dynamics friction coefficient(correspond to gamma_ln in AMBER)
ld_seed = 1732  ; random seed

; Neighborsearching
ns_type         = grid          ; search neighboring grid cells
nstlist         = 20              ; 10 fs
; rlist           = 1.0 ; short-range neighborlist cutoff (in nm)  ; 1 is the default
; rcoulomb        = 1.0 ; short-range electrostatic cutoff (in nm) ; 1 is the default
; rvdw            = 1.0 ; short-range van der Waals cutoff (in nm) ; 1 is the default

; Electrostatics
coulombtype     = PME    ; Particle Mesh Ewald for long-range electrostatics
; pme_order       = 4    ; cubic interpolation : 4 is the default
; fourierspacing  = 0.12 ; grid spacing for FFT(nm) ; 0.12 is the default

; Temperature coupling is on
tcoupl    = V-rescale   ; modified Berendsen thermostat
tc-grps   = Non-Water Water   ; two coupling groups - more accurate
tau_t     = 0.1  0.1               ; time constant, in ps
ref_t     = 300  300      ; reference temperature, one for each group, in K

; Pressure coupling is on
pcoupl         = Parrinello-Rahman      ; Pressure coupling on in NPT
pcoupltype     = isotropic              ; uniform scaling of box vectors
tau_p          = 2.0              ; time constant, in ps
ref_p          = 1.0                ; reference pressure, in bar
compressibility = 4.5e-5  ; isothermal compressibility of water, bar^-1
refcoord_scaling = com

; Periodic boundary conditions
pbc                 = xyz               ; 3-D PBC

; Dispersion correction
DispCorr     = EnerPres ; account for cut-off vdW scheme

; Velocity generation
gen_vel         = no            ; Velocity generation
gen_temp        = 300.0         ; temperature for Maxwell distribution
gen_seed        = 1732          ; used to initialize random generator for random velocities
"""


def make_gromacs_directories(parent_dir):
    pathlib.Path("%s/top" % parent_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path("%s/minimize" % parent_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path("%s/heat" % parent_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path("%s/pr" % parent_dir).mkdir(parents=True, exist_ok=True)
    return parent_dir, "%s/top" % parent_dir, "%s/minimize" % parent_dir, \
        "%s/heat" % parent_dir, "%s/pr" % parent_dir


def copy_gromacs_files(inputdir, topdir, name):
    shutil.copy("%s/%s.gro" % (inputdir, name), "%s/%s.gro" % (topdir, name))
    shutil.copy("%s/%s.top" % (inputdir, name), "%s/%s.top" % (topdir, name))
    shutil.copy("%s/index.ndx" % inputdir, "%s/index.ndx" % topdir)

    subprocess.getoutput("cp %s/*.top %s/*.itp %s" % (inputdir, inputdir, topdir))


def gen_minimize_mdp(defines, filename, steps=200, emstep=0.01):
    if not defines:
        defines = ""

    with open(filename, "w") as fout:
        fout.write(template_minimize.format(defines=defines, steps=steps, emstep=emstep))


def gen_heat_mdp(nsteps, defines, path,
                 dt=0.002, pcoupl="berendsen",
                 gen_vel=False, annealing=False,
                 continuation="yes"):
    if not defines:
        defines = ""

    if gen_vel:
        velocity_setting = "gen_vel = yes\ngen_seed = -1"
    else:
        velocity_setting = "gen_vel = no"
    if annealing:
        annealing_setting = """
annealing = single single
annealing_npoints = 2 2
annealing_time = 0 200 0 200
annealing_temp = 0 300 0 300
"""
    else:
        annealing_setting = ""

    with open(path, "w") as fout:
        fout.write(template_heat % (defines, dt, nsteps, pcoupl,
                                    velocity_setting, continuation,
                                    annealing_setting))


def gen_pr_mdp(nsteps, defines, snapshot_interval, path):
    with open(path, "w") as fout:
        fout.write(template_pr.format(nsteps=nsteps,
                                      defines=defines,
                                      snapshot_interval=snapshot_interval))


def gen_mdrun_job(template_file, name, path, post_comm=""):
    data = {
        "NAME": name,
        "POST_COMMAND": post_comm
    }
    env = jinja2.Environment(loader=jinja2.FileSystemLoader("/"))
    template = env.get_template(template_file)
    with open(path, "w") as fout:
        fout.write(template.render(data))


def update_config(data, var_dict):
    # TODO check key validity
    for key, value in var_dict.items():
        key1, key2 = key.split(":")
        if key1 not in data.keys():
            data[key1] = {}
        data[key1][key2] = value


def validate_config(data):
    # TODO check more things
    if("calc_dir" in data["General"]) and ("output_dir" not in data["General"]):
        print("WARNING: variable General:calc_dir is deprecated. Use General:output_dir instead.")
        data["General"]["output_dir"] = data["General"]["calc_dir"]
    assert("output_dir" in data["General"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run gromacs jobs automatically")
    parser.add_argument("--version", action="version", version=VERSION)

    # https://gist.github.com/vadimkantorov/37518ff88808af840884355c845049ea
    parser.add_argument("-v", default={},
                        action=type('', (argparse.Action, ),
                                    dict(__call__=lambda a, p, n, v, o:
                                         getattr(n, a.dest).update(dict([v.split('=')])))))

    parser.add_argument("conf", help="configfile")
    parser.add_argument("tmpl_mdrun")
    args = parser.parse_args()

    dat = configparser.ConfigParser()
    dat.read_dict({
        "ProductionRun": {
            "snapshot_interval": "5000",
            "defines": "",
            "post_comm": "",
        },
    })
    dat.read(args.conf, "UTF-8")
    update_config(dat, args.v)
    validate_config(dat)
    # print(dat["General"]["name"])

    # 1. make directories
    PARENT_DIR, TOP_DIR, MIN_DIR, HEAT_DIR, PR_DIR\
        = make_gromacs_directories(dat["General"]["output_dir"])

    # 1.1 copy raw data
    copy_gromacs_files(dat["General"]["input_dir"], TOP_DIR,
                       dat["General"]["name"])

    # 2. MINIMIZE
    # min1.mdp : minimization for hydrogen atoms
    # min2.mdp : minimization for all atoms
    gen_minimize_mdp("-DPOSRES1000", "%s/min1.mdp" % MIN_DIR)
    gen_minimize_mdp("",             "%s/min2.mdp" % MIN_DIR)

    # 3. HEAT
    gen_heat_mdp(nsteps=100000, defines="-DPOSRES1000",
                 gen_vel=True, annealing=True, pcoupl="no",
                 continuation="no",
                 path="%s/md1.mdp" % HEAT_DIR)
    gen_heat_mdp(nsteps=50000, defines="-DPOSRES1000", path="%s/md2.mdp" % HEAT_DIR)
    gen_heat_mdp(nsteps=50000, defines="-DPOSRES500",  path="%s/md3.mdp" % HEAT_DIR)
    gen_heat_mdp(nsteps=50000, defines="-DPOSRES200",  path="%s/md4.mdp" % HEAT_DIR)
    gen_heat_mdp(nsteps=50000, defines="-DPOSRES100",  path="%s/md5.mdp" % HEAT_DIR)
    gen_heat_mdp(nsteps=50000, defines="-DPOSRES50",   path="%s/md6.mdp" % HEAT_DIR)
    gen_heat_mdp(nsteps=50000, defines="-DPOSRES20",   path="%s/md7.mdp" % HEAT_DIR)
    gen_heat_mdp(nsteps=50000, defines="-DPOSRES10",   path="%s/md8.mdp" % HEAT_DIR)
    gen_heat_mdp(nsteps=50000, defines="-DPOSRES0",    path="%s/md9.mdp" % HEAT_DIR)

    # 4. PRODUCTION RUN
    gen_pr_mdp(nsteps=dat["ProductionRun"]["steps_pr"],
               defines=dat["ProductionRun"]["defines"],
               snapshot_interval=dat["ProductionRun"]["snapshot_interval"],
               path="{}/{}.mdp".format(PR_DIR, dat["General"]["name"]))

    gen_mdrun_job(args.tmpl_mdrun,
                  dat["General"]["name"],
                  "%s/mdrun.sh" % PARENT_DIR,
                  dat["ProductionRun"]["post_comm"])
