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

finished_info=finished_step_list
touch $finished_info

for stepname in {{ STEP_NAMES }}
do
    prev=$now
    now=$stepname

    if [ `grep -x ${now} $finished_info | wc -l` = 1 ] ;then
       continue
    fi
    
    rm -f ${now}.out.mdp ${now}.tpr ${now}.log ${now}.gro ${now}.trr ${now}.edr ${now}.cpt
    echo $GMX grompp -f ${now}.mdp -o ${now}.tpr \
      -c ${prev}.gro -p ${top} \
      -r ${prev}.gro -n index.ndx
    echo $GMX mdrun -reprod -v -s ${now}.tpr \
      -cpo ${now}.cpt -x ${now}.xtc -c ${now}.gro -e ${now}.edr -g ${now}.log \
      && echo $now >> $finished_info

    $GMX grompp -f ${now}.mdp -o ${now}.tpr \
      -c ${prev}.gro -p ${top} \
      -r ${prev}.gro -n index.ndx
    $GMX mdrun -reprod -v -s ${now}.tpr \
      -cpo ${now}.cpt -x ${now}.xtc -c ${now}.gro -e ${now}.edr -g ${now}.log \
      && echo $now >> $finished_info

done


ln -s $stepname.xtc {{ OUT_TRAJ }}
