logging(){
    echo [`date`] $1
}
logging_error(){
    logging "[ERROR] $1"
}
logging_info(){
    logging "[INFO]  $1"
}
logging_debug(){
    logging "[DEBUG] $1"
}
logging_warn(){
    logging "[WARN]  $1"
}

is_calculated(){
    local CALCDIR=$1
    local i=$2
    local mapfile=$3
    if [ -f $CALCDIR/system$i/$mapfile ]; then
	echo 1
    else
	echo 0
    fi
}

get_gpu_list(){
    local gpu_lst=(${CUDA_VISIBLE_DEVICES//,/ })
    echo ${gpu_lst[@]}
}

get_ncpus(){
    grep processor /proc/cpuinfo | wc -l
}

#https://qiita.com/srea/items/28073bc90d65eed0856d
get_ini_variable(){
    local ini_file=$1
    local ini_section=$2
    local ini_variable=$3

    if [ ! -f $ini_file ];then
	echo "ini file does not found" 1>&2
	exit 1
    fi

    sed -e 's/[[:space:]]*\=[[:space:]]*/=/g' \
	-e 's/;.*$//' \
	-e 's/[[:space:]]*$//' \
	-e 's/^[[:space:]]*//' \
	-e "s/^\(.*\)=\([^\"']*\)$/\1=\2/" \
	< $ini_file \
    	| sed -n -e "/^\[$ini_section\]/,/^\s*\[/{/^[^;].*\=.*/p;}" \
	| grep ^$ini_variable | cut -d= -f2
}
