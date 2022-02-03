iter=20
iter_st=0
INPUTDIR=`pwd`
OUTPUTDIR=`pwd`/output/test
map_prefix=PMAP
protein_param_file=example/protein.conf
probe_param_file=example/probe.conf
TARGET_NAME=TEST_PROJECT
steps_pr=20000000          # 40 ns (timestep: 2fs)
snapshot_interval=5000     # 10 ps per snapshot
begin_snapshot_for_pr=2001 # 20--40 ns

# NCPUS=14
# CUDA_VISIBLE_DEVICES="0"
