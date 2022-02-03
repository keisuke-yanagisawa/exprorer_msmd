#! /usr/bin/python3

import argparse
import configparser
import pathlib
import subprocess
import shutil
import jinja2
import yaml

VERSION = "1.0.0"

template_minimize = """
define       = {define}
;
integrator   = steep
nsteps       = {nsteps}
nstlog       = {nstlog}
;
nstlist      = 5
ns_type      = grid
pbc          = {pbc}
;
coulombtype  = PME
;
vdwtype      = cut-off
;
constraints  = none
"""

template_heat = """
define       = {define}
;
integrator   = md
dt           = {dt}
nsteps       = {nsteps}
;
nstxtcout     = {nstxtcout}
xtc-precision = 10000
nstlog        = {nstlog}
nstenergy     = {nstenergy}
;
nstlist      = 5
ns_type      = grid
pbc          = {pbc}
;
coulombtype  = PME
;
vdwtype      = cut-off
dispcorr     = EnerPres
;
tcoupl       = v-rescale
tc_grps      = !Water_and_ions Water_and_ions
tau_t        = 0.1 0.1
ref_t        = {target_temp} {target_temp}
;
pcoupl       = {pcoupl}
pcoupltype   = isotropic
compressibility = 4.5e-5
ref_p        = {pressure}
refcoord_scaling = all
;
gen_vel = yes
gen_seed = {seed}
gen_temp = {initial_temp}
;
continuation = no
constraints  = hbonds
constraint_algorithm = LINCS
;
Annealing = single single
Annealing_npoints = 2 2
Annealing_time = 0 {nsteps*dt} 0 {nsteps*dt}
Annealing_temp = {initial_temp} {target_temp} {initial_temp} {target_temp}
"""

template_equil = """
define       = {define}
;
integrator   = md
dt           = {dt}
nsteps       = {nsteps}
;
nstxtcout     = {nstxtcout}
xtc-precision = 10000
nstlog        = {nstlog}
nstenergy     = {nstenergy}
;
nstlist      = 5
ns_type      = grid
pbc          = {pbc}
;
coulombtype  = PME
;
vdwtype      = cut-off
dispcorr     = EnerPres
;
tcoupl       = v-rescale
tc_grps      = !Water_and_ions Water_and_ions
tau_t        = 0.1 0.1
ref_t        = {temperature} {temperature}
;
pcoupl       = {pcoupl}
pcoupltype   = isotropic
compressibility = 4.5e-5
ref_p        = {pressure}
refcoord_scaling = all
;
gen_vel = no
;
continuation = yes
constraints  = hbonds
constraint_algorithm = LINCS
"""

template_pr = """
define            = {define}
;
integrator      = md
dt              = {dt}
nsteps          = {nsteps}

; Output control
nstxtcout       = {nstxtcout}          ; save coordinates every N ps
xtc-precision   = 10000
nstenergy       = {nstenergy}
nstlog          = {nstlog}

; Bond parameters
continuation    = yes           ; Restarting after previous RUN
constraint_algorithm = LINCS    ; holonomic constraints
constraints          = hbonds  ; Convert the bonds with H-atoms to constraints
; Langevin Dynamics
bd_fric = 2.0           ; Brownian dynamics friction coefficient(correspond to gamma_ln in AMBER)
ld_seed = {seed}
;
ns_type         = grid
nstlist         = 5
;
coulombtype     = PME    ; Particle Mesh Ewald for long-range electrostatics
;
tcoupl    = V-rescale
tc-grps   = !Water_and_ions Water_and_ions
tau_t     = 0.1  0.1
ref_t     = {temperature} {temperature}
; Pressure coupling is on
pcoupl         = {pcoupl}
pcoupltype     = isotropic
ref_p          = {pressure}
compressibility = 4.5e-5
refcoord_scaling = com
; Periodic boundary conditions
pbc                 = {pbc}
; Dispersion correction
DispCorr     = EnerPres ; account for cut-off vdW scheme
; Velocity generation
gen_vel         = no            ; Velocity generation
"""


def make_gromacs_directories(parent_dir):
    pathlib.Path("%s/top" % parent_dir).mkdir(parents=True, exist_ok=True)
    pathlib.Path("%s/simulation" % parent_dir).mkdir(parents=True, exist_ok=True)
    return parent_dir, "%s/top" % parent_dir, "%s/simulation" % parent_dir


def copy_gromacs_files(inputdir, topdir, name):
    shutil.copy("%s/%s.gro" % (inputdir, name), "%s/%s.gro" % (topdir, name))
    shutil.copy("%s/%s.top" % (inputdir, name), "%s/%s.top" % (topdir, name))
    shutil.copy("%s/index.ndx" % inputdir, "%s/index.ndx" % topdir)

    # TODO#19 detect copy files automatically (from top file)
    for i in [0, 10, 20, 50, 100, 200, 500, 1000]:
        shutil.copy("%s/posrePROTEIN%d.itp" % (inputdir, i),
                    "%s/posrePROTEIN%d.itp" % (topdir, i))
    subprocess.getoutput("cp %s/*.top %s/*.itp %s" % (inputdir, inputdir, topdir))

def gen_mdp(protocol_dict, MD_DIR):
    if protocol_dict["step"] == "minimization":
        template = template_minimize
    elif protocol_dict["step"] == "heating":
        template = template_heat
        if "target_temp" not in protocol_dict:
            protocol_dict["target_temp"] = protocol_dict["temperature"]
        if "initial_temp" not in protocol_dict:
            protocol_dict["initial_temp"] = 0
        protocol_dict["nsteps*dt"] = protocol_dict["nsteps"] * protocol_dict["dt"]
    elif protocol_dict["step"] == "equilibration":
        template = template_equil
    elif protocol_dict["step"] == "production":
        template = template_pr

    with open(f"{MD_DIR}/{protocol_dict['name']}.mdp", "w") as fout:
        fout.write(template.format(**protocol_dict))

def gen_mdrun_job(template_file, step_names, name, path, post_comm=""):
    data = {
        "NAME"         : name,
        "POST_COMMAND" : post_comm,
        "STEP_NAMES"   : " ".join(step_names)
    }

    env = jinja2.Environment(loader=jinja2.FileSystemLoader("."))
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

### protocol_dict:define must be rewritten to "" if it is None
    

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
    dat.read(args.conf, "UTF-8")
    update_config(dat, args.v)
    validate_config(dat)
    # print(dat["General"]["name"])

    with open("/gs/hs0/tga-pharma/yanagisawa/workspace/9999_git_repositories/exprorer/cmd_calculation/example/md_protocol.yaml") as fin:
        protocol = yaml.safe_load(fin)

    # 1. make directories
    PARENT_DIR, TOP_DIR, MD_DIR = make_gromacs_directories(dat["General"]["output_dir"])

    for i in range(len(protocol["sequence"])):
        step = protocol["sequence"][i]
        if "name" not in step:
            step["name"] = f"step{i+1}"

        if step["define"] is None:
            step["define"] = ""
        
        step.update(protocol["general"])
        gen_mdp(step, MD_DIR)
        protocol["sequence"][i] = step # update
        

    gen_mdrun_job(args.tmpl_mdrun,
                  [d["name"] for d in protocol["sequence"]],
                  dat["General"]["name"],
                  "%s/mdrun.sh" % PARENT_DIR)#,
#                  dat["ProductionRun"]["post_comm"])

    # 1.1 copy raw data
    copy_gromacs_files(dat["General"]["input_dir"], TOP_DIR,
                       dat["General"]["name"])

