#! /bin/bash

#### init ####

. setting/initialize.sh

#### utility functions ####

. script/bash_functions.sh

#### calculation functions ####

genpmap(){
    local i=$1
    $PYTHON $WORKDIR/script/genpmap_main.py -basedir $ANALYSISDIR/system$i \
        $protein_param_file $probe_param_file --debug
}

maxpmap(){
    #TODO: map_test_**.dx_pmap.dx is terrible hard-coding
    $PYTHON $WORKDIR/script/maxpmap_main.py $ANALYSISDIR/system*/map_test_nV.dx_pmap.dx \
        $ANALYSISDIR/maxPMAP_nV.dx
    $PYTHON $WORKDIR/script/maxpmap_main.py $ANALYSISDIR/system*/map_test_nVH.dx_pmap.dx \
        $ANALYSISDIR/maxPMAP_nVH.dx
    $PYTHON $WORKDIR/script/maxpmap_main.py $ANALYSISDIR/system*/map_test_O.dx_pmap.dx \
        $ANALYSISDIR/maxPMAP_V.dx
    $PYTHON $WORKDIR/script/maxpmap_main.py $ANALYSISDIR/system*/map_test_V.dx_pmap.dx \
        $ANALYSISDIR/maxPMAP_O.dx
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

for i in `seq $iter_st $(( $iter_ed - 1 ))`
do
    genpmap $i &
done
wait

maxpmap
