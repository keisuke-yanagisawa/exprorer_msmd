#!/bin/sh

ncpus=$1

## initialize
top={{ TOP }}
now={{ GRO | replace(".gro", "") }}
if [ A$GMX = "A" ];then
    GMX=gmx
fi

## for OpenMP parallelization
export OMP_NUM_THREADS=$ncpus
## for thread-MPI parallelization

for stepname in {{ STEP_NAMES }}
do
    prev_cpt=""
    prev=$now
    now=$stepname

    rm -f ${now}.out.mdp ${now}.tpr ${now}.log ${now}.gro ${now}.trr ${now}.edr ${now}.cpt
    echo $GMX grompp -f ${now}.mdp -po ${now}.out.mdp -o ${now}.tpr \
	-c ${prev}.gro -p ${top} \
	-r ${prev}.gro -n index.ndx ${prev_cpt}
    echo $GMX mdrun -reprod -v -s ${now} -deffnm ${now} -cpi ${now}.cpt

    $GMX grompp -f ${now}.mdp -po ${now}.out.mdp -o ${now}.tpr \
	  -c ${prev}.gro -p ${top} \
	  -r ${prev}.gro -n index.ndx ${prev_cpt}
    $GMX mdrun -reprod -v -s ${now} -deffnm ${now} -cpi ${now}.cpt

    prev_cpt="-t $now.cpt"
done

# trr -> cpt
# -t ${prev_cpt} trr? cpt?

ln -s $stepname.xtc {{ OUT_TRAJ }}