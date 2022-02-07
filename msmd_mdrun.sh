#! /bin/bash

#### init ####

. setting/initialize.sh

if [ $do_reorder"A" == "A" ] ; then
    do_reorder="Y"
fi

#### utility functions ####

. script/bash_functions.sh

#### calculation functions ####

preparation(){

    local i=$1
    local do_reorder=$2 # Y or otherwise (N)
    cd $OUTPUTDIR/prep$i

    if [ `is_calculated $OUTPUTDIR $i ${map_prefix}_nVH.dx` == 1 ]; then
        logging_info "$TARGET_NAME $i : skipped because it has been already calculated"
        continue
    fi

    $GMX make_ndx -f $OUTPUTDIR/prep$i/$TARGET_NAME.gro << EOF
q
EOF

    $PYTHON $WORKDIR/script/mdrun.py \
	DUMMY.conf \
	$WORKDIR/script/mdrun.sh \
	-v General:input_dir=$OUTPUTDIR/prep$i \
	-v General:output_dir=$OUTPUTDIR/system$i \
	-v General:name=$TARGET_NAME \
    -v General:protocol_yaml=$OUTPUTDIR/prep$i/input/md_protocol.yaml \
	-v ProductionRun:steps_pr=$steps_pr \
	-v ProductionRun:snapshot_interval=$snapshot_interval
}

start_mdrun(){
    local TARGET_NAME=$1
    local i=$2
    local ncpus=$3
    
    if [ `is_calculated $OUTPUTDIR $i ${map_prefix}_nVH.dx` == 1 ]; then
        logging_info "$TARGET_NAME $i : skipped because it has been already calculated"
        return 0
    fi

    cd $OUTPUTDIR/system$i
    GMX=$GMX bash mdrun.sh $ncpus
}


#### output environments ####

logging_info "which $PYTHON: `which $PYTHON`"

#### main routine ####

if [ $# != 1 ]; then
    logging_debug "$# != 1"
    logging_error "Invalid usage has been detected!"
    logging_info  "Usage: bash $0 CONFIG.sh"
    exit 1
fi

CONFIG_FILE=$1
source $CONFIG_FILE

#### preparation ####

WORKDIR=`pwd`
unset OMP_NUM_THREADS
NCPUS=${NCPUS:-`get_ncpus`}

gpus=( `get_gpu_list` )
echo $gpus
if [ ${#gpus[@]} == "0" ]; then
    gpus=0 # dummy information
    # logging_error "There is no available GPU according to the CUDA_VISIBLE_DEVICES variable."
    # logging_error "This code assume at least 1 GPU. Please set the variable correctly."
    # logging_error "Exit."
    # exit 1
fi

NCPUS_PER_GPU=$(( $NCPUS / ${#gpus[@]} ))

iter_ed=$(( $iter + $iter_st ))

# parallel preparation
for i in `seq $iter_st $(( $iter_ed - 1 ))`
do
    preparation $i $do_reorder &
done
wait

echo parallel execution
# parallel execution
i=$iter_st
while [ $i -lt $iter_ed ]
do
    for j in ${gpus[@]}
    do
	echo $i $iter_st $iter_ed "(" ${gpus[@]} ")" $j $NCPUS_PER_GPU
	export CUDA_VISIBLE_DEVICES="$j"
	start_mdrun $TARGET_NAME $i $NCPUS_PER_GPU &
	
	i=$(( i + 1 ))
	if [ $i -ge $iter_ed ]; then
	    break
	fi
    done
    wait
done

