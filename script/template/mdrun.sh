#!/bin/sh

hostname
ncpus=$1

## initialize
top={{ TOP }}
now={{ GRO | replace(".gro", "") }}
if [ A$GMX = "A" ];then
    GMX=gmx
fi


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
    echo $GMX grompp -maxwarn 1 -f ${now}.mdp -o ${now}.tpr \
      -c ${prev}.gro -p ${top} \
      -r ${prev}.gro -n index.ndx
    echo $GMX mdrun -nt $ncpus -v -s ${now}.tpr \
      -cpo ${now}.cpt -x ${now}.xtc -c ${now}.gro -e ${now}.edr -g ${now}.log

    $GMX grompp -maxwarn 1 -f ${now}.mdp -o ${now}.tpr \
      -c ${prev}.gro -p ${top} \
      -r ${prev}.gro -n index.ndx
    $GMX mdrun -nt $ncpus -v -s ${now}.tpr \
      -cpo ${now}.cpt -x ${now}.xtc -c ${now}.gro -e ${now}.edr -g ${now}.log \
      || exit

    echo $now >> $finished_info

done


ln -s $stepname.xtc {{ OUT_TRAJ }}
