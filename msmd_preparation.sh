#! /bin/bash

#### init ####

. setting/initialize.sh

#### utility functions ####

. script/bash_functions.sh

#### calculation functions ####

preparation(){

    local i=$1

    if [ `is_calculated $OUTPUTDIR $i ${map_prefix}_nVH.dx` == 1 ]; then
        logging_info "$TARGET_NAME $i : skipped because it has been already calculated"
        continue
    fi

    rm -r $OUTPUTDIR/prep$i # remove existing files
    mkdir -p $OUTPUTDIR/prep$i/input

    logging_debug "protein_param_file $protein_param_file"
    logging_debug "probe_param_file $probe_param_file"

    cp $protein_param_file $OUTPUTDIR/prep$i/input/protein.conf
    cp $probe_param_file   $OUTPUTDIR/prep$i/input/probe.conf
    cd $OUTPUTDIR/prep$i

    cosolvent_ID=`get_ini_variable $OUTPUTDIR/prep$i/input/probe.conf Cosolvent cid`
    $PYTHON $WORKDIR/script/cosolvent_box_generation.py \
	-prot_param $OUTPUTDIR/prep$i/input/protein.conf \
	-cosolv_param $OUTPUTDIR/prep$i/input/probe.conf \
	-tin   $WORKDIR/script/template_leap.in \
	-oprefix $OUTPUTDIR/prep$i/$TARGET_NAME \
        --packmol $PACKMOL \
	--tleap   $TLEAP \
	-seed $i \
	-no-rm-temp

    $PYTHON $WORKDIR/script/convert_amber_to_gromacs.py \
	   -iprefix $OUTPUTDIR/prep$i/${TARGET_NAME} \
	   -oprefix $OUTPUTDIR/prep$i/${TARGET_NAME}_GMX

    # add virtual interaction between cosolvents
    $PYTHON $WORKDIR/script/addvirtatom2top.py \
	-i $OUTPUTDIR/prep$i/${TARGET_NAME}_GMX.top \
	-o $OUTPUTDIR/prep$i/${TARGET_NAME}_tmp.top \
	-cname ${cosolvent_ID} \
	-ovis $OUTPUTDIR/prep$i/virtual_repulsion.top 
    $PYTHON $WORKDIR/script/addvirtatom2gro.py \
	-i $OUTPUTDIR/prep$i/${TARGET_NAME}_GMX.gro \
	-o $OUTPUTDIR/prep$i/$TARGET_NAME.gro \
	-cname ${cosolvent_ID} \
	-vname "VIS"

    # gen position restraint files
    $PYTHON $WORKDIR/script/add_posredefine2top.py \
	-v -res WAT Na+ Cl- CA MG ZN CU ${cosolvent_ID} \
	-target protein \
	-gro $OUTPUTDIR/prep$i/$TARGET_NAME.gro \
	-i $OUTPUTDIR/prep$i/${TARGET_NAME}_tmp.top \
	-o $OUTPUTDIR/prep$i/$TARGET_NAME.top \


    $GMX trjconv -s $OUTPUTDIR/prep$i/$TARGET_NAME.gro -f $OUTPUTDIR/prep$i/$TARGET_NAME.gro -o $OUTPUTDIR/prep$i/$TARGET_NAME.pdb <<EOF
0
EOF
    # prepared files are:
    # xxxx.top, xxxx.gro, virtual_repulsion.top, xxxx.pdb
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

iter_ed=$(( $iter + $iter_st ))

# parallel preparation
for i in `seq $iter_st $(( $iter_ed - 1 ))`
do
    preparation $i &
done
wait
