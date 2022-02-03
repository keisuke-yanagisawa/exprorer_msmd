. $HOME/.tsubame3.sh

module load amber/18up17
export PATH=$(echo -n $PATH |tr ':' '\n' |sed "/\/apps\/t3\/sles12sp2\/isv\/amber\/amber18up17\/anaconda3\/bin/d" |tr '\n' ':')
module load cuda/9.2.148
module load gromacs/2019.4

export CUDA_VISIBLE_DEVICES=0,1,2,3
export PYTHONPATH="$WORKSPACE/9999_git_repositories/ak-library/python:$PYTHONPATH"

conda activate md_analysis

PACKMOL=packmol
TLEAP=tleap
CPPTRAJ=cpptraj
GMX=gmx_mpi
PYTHON=python
