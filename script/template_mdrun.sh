#!/bin/sh

ncpus=$1

## initialize
top_wo_ext=../top/{{ NAME }}
now=../top/{{ NAME }}
if [ A$GMX = "A" ];then
    GMX=gmx
fi

## for OpenMP parallelization
export OMP_NUM_THREADS=$ncpus
## for thread-MPI parallelization
#tMPI_COMM="-nt $ncpus"

## minimize
cd minimize

for i in 1 2
do
    prev=$now
    now=../minimize/min$i

    rm -f ${now}.out.mdp ${now}.tpr ${now}.log ${now}.gro ${now}.trr ${now}.edr ${now}.cpt
    $GMX grompp -f ${now}.mdp -po ${now}.out.mdp -o ${now}.tpr \
	-c ${prev}.gro -p ${top_wo_ext}.top \
	-r ${prev}.gro
    $GMX mdrun -v -s ${now} -deffnm ${now} $tMPI_COMM     
done
cd ..

## heating - equillibration
cd heat
for i in 1 2 3 4 5 6 7 8 9
do
    prev=$now
    if [ $i = 1 ];then
	prev_cpt=$now.trr
    else
	prev_cpt=$now.cpt
    fi
    now=../heat/md$i

    rm -f ${now}.out.mdp ${now}.tpr ${now}.log ${now}.gro ${now}.trr ${now}.edr ${now}.cpt
    $GMX grompp -f ${now}.mdp -po ${now}.out.mdp -o ${now}.tpr \
	-c ${prev}.gro -p ${top_wo_ext}.top \
	-r ${prev}.gro -n ../top/index.ndx -t ${prev_cpt}
    $GMX mdrun -v -s ${now} -deffnm ${now} -cpi ${now}.cpt $tMPI_COMM 
done
cd ..

## production run
cd pr
prev=$now
now=../pr/{{ NAME }}
rm -f ${now}.out.mdp ${now}.tpr ${now}.log ${now}.gro ${now}.trr ${now}.edr ${now}.cpt
$GMX grompp -f ${now}.mdp -po ${now}.out.mdp -o ${now}.tpr \
    -c ${prev}.gro -p ${top_wo_ext}.top \
    -r ${prev}.gro -n ../top/index.ndx -t ${prev}.cpt
$GMX mdrun -v -s ${now} -deffnm ${now} -cpi ${now}.cpt $tMPI_COMM 
cd ..

## post calculation
cd pr
{{ POST_COMMAND }}
cd ..
