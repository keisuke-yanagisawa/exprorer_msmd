import argparse
import os
from pathlib import Path
from string import Template

import yaml

from script.utilities import util
from script.utilities.logger import logger


def read_yaml(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MSMD simulation")
    parser.add_argument("setting_yaml")
    parser.add_argument("-v,--verbose", dest="verbose", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--skip-preprocess", action="store_true")
    parser.add_argument("--skip-simulation", action="store_true")
    parser.add_argument("--skip-postprocess", action="store_true")
    parser.add_argument("-N,--run-per-job", dest="run_per_job", type=int, default=1)
    args = parser.parse_args()

    environment = read_yaml(os.path.dirname(util.getabsolutepath(Path(__file__))) + "/environment.yaml")

    # initial logger level is "warn"
    if args.debug:
        logger.setLevel("debug")
    elif args.verbose:
        logger.setLevel("info")

    logger.info(f"read yaml: {args.setting_yaml}")
    setting = util.parse_yaml(args.setting_yaml)
    indices = list(set(util.expand_index(setting["general"]["iter_index"])))
    workdir = setting["general"]["workdir"]

    GROUP = ""
    QTYPE = ""
    MAXTIME = ""
    if environment["use_scheduler"] == True:
        GROUP = environment["scheduler"]["group"]
        QTYPE = environment["scheduler"]["qtype"]
        MAXTIME = environment["scheduler"]["maxtime"]

    use_singularity = environment["use_singularity"]
    singularity_sifpath = ""
    singularity_bind = ""
    singularity_prerequirement = ""
    if use_singularity == True:
        singularity_sifpath = environment["singularity"]["sifpath"]
        singularity_prerequirement = environment["singularity"]["prerequirement"]
        singularity_bind = environment["singularity"]["bind"]

    template = ""
    if environment["use_scheduler"] != True:
        template = open(os.path.dirname(util.getabsolutepath(Path(__file__))) + "/template/pbs", "r").read()
    elif environment["scheduler"]["type"] == "PBS":
        template = open(os.path.dirname(util.getabsolutepath(Path(__file__))) + "/template/pbs", "r").read()
    elif environment["scheduler"]["type"] == "Slurm":
        template = open(os.path.dirname(util.getabsolutepath(Path(__file__))) + "/template/slurm", "r").read()

    import math

    njobs = math.ceil(len(indices) / args.run_per_job)
    for i, gr in enumerate([indices[i::njobs] for i in range(njobs)]):
        JOB_NAME = setting["general"]["name"]
        ITER_INDEX = ",".join([str(j) for j in gr])
        SETTING_YAML = args.setting_yaml

        # generate slurm script
        os.system(f"mkdir -p {workdir}")
        with open(f"{workdir}/JOB{i}.sh", "w") as fout:
            fout.write(
                template.format(
                    **{
                        "GROUP": GROUP,
                        "QTYPE": QTYPE,
                        "runID": f"{i}",
                        "JOB_NAME": JOB_NAME,
                        "TIME_LENGTH": MAXTIME,
                        "use_singularity": use_singularity,
                        "singularity_prerequirement": singularity_prerequirement,
                        "singularity_sifpath": singularity_sifpath,
                        "singularity_bind": singularity_bind,
                        "PATH_EXPRORER_MSMD": os.path.dirname(util.getabsolutepath(Path(__file__))) + "/exprorer_msmd",
                        "ITER_INDEX": ITER_INDEX,
                        "YAML": SETTING_YAML,
                        "NGPUS": len(gr),
                    }
                )
            )
        os.system(f"sbatch {workdir}/JOB{i}.sh")
