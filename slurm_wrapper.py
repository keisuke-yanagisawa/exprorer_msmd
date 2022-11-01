import argparse
import os
from script.utilities.logger import logger
from script.utilities import util

slurm_template="""#!/usr/bin/bash
#SBATCH -J {JOB_NAME}
#SBATCH -p {PARTITION}
#SBATCH -N 1
#SBATCH --gres=gpu:1
#SBATCH --get-user-env
#SBATCH -o {JOB_NAME}.out
#SBATCH -e {JOB_NAME}.err

EXPRORER_MSMD="/mnt/fs/yanagisawa/workspace/9999_git_repositories/exprorer_msmd/exprorer_msmd"
MSMD_OPTIONS="--iter-index {ITER_INDEX}"
SETTING_YAML={SETTING_YAML}

#IF SINGULARITY
module load singularity
SINGULARITY_OPTIONS="--bind /mnt/ghosts:/mnt/ghosts --bind /mnt/fs:/mnt/fs --bind /tmp:/tmp"
SINGULARITY_SIF="/mnt/fs/yanagisawa/workspace/9999_git_repositories/exprorer_msmd/.devcontainer/Singularity.sif"
SINGULARITY="singularity exec --nv $SINGULARITY_OPTIONS $SINGULARITY_SIF"

echo $SINGULARITY $EXPRORER_MSMD $MSMD_OPTIONS $SETTING_YAML
$SINGULARITY $EXPRORER_MSMD $MSMD_OPTIONS $SETTING_YAML
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run MSMD simulation")
    parser.add_argument("setting_yaml")
    parser.add_argument("-v,--verbose", dest="verbose", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--skip-preprocess", action="store_true")
    parser.add_argument("--skip-simulation", action="store_true")
    parser.add_argument("--skip-postprocess", action="store_true")
    parser.add_argument("-N,--run-per-job", dest="run_per_job", type=int, default=1)
    args = parser.parse_args()

    # initial logger level is "warn"
    if args.debug:
        logger.setLevel("debug")
    elif args.verbose:
        logger.setLevel("info")

    logger.info(f"read yaml: {args.setting_yaml}")
    setting = util.parse_yaml(args.setting_yaml)
    indices = list(set(util.expand_index(setting["general"]["iter_index"])))    
    workdir = setting["general"]["workdir"]

    import math
    njobs = math.ceil(len(indices) / args.run_per_job)
    for i, gr in enumerate([indices[i::njobs] for i in range(njobs)]):
        JOB_NAME = setting["general"]["name"]
        PARTITION = "defq"
        ITER_INDEX = ",".join([str(i) for i in gr])
        SETTING_YAML = args.setting_yaml

        # generate slurm script
        with open(f"{workdir}/slurm{i}.sh", "w") as fout:
            fout.write(slurm_template.format(**{
                "JOB_NAME": f"JOB{i}_{JOB_NAME}",
                "PARTITION": PARTITION,
                "ITER_INDEX": ITER_INDEX,
                "SETTING_YAML": SETTING_YAML,
            }))
        os.system(f"sbatch {workdir}/slurm{i}.sh")
