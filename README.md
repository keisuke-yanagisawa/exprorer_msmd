# EXPRORER_MSMD

MSMD (mixed-solvent molecular dynamics) engine and analysis tools

## Requirements

Below requirements are based on our environment

- python 3 (`$PYTHON`)
  - Several modules are needed (an example of installation is shown below).
- openbabel 2.4.1
  - **Do not use the latest version, 3.0.0 or newer** since we have not tested yet.
  - It is possible to install it along with python 3 via conda.
- AmberTools 18 (`$TLEAP`, `$CPPTRAJ`)
  - It is needed for use of `tleap` and `cpptraj`
- Gromacs 2019.4 (`$GMX`)
- packmol 18.169 (`$PACKMOL`)
- (Gaussian 16.B01)
  - So far we have not included the codes for QM calculation of probes (structural optimization and RESP charge calculation)
  - It is needed **only** when new probes are introduced.

## Preparation

It can be used with Docker. Below instruction is for user who will not use Docker.

### Construction of Python3 environment with Conda

```sh
conda create -n ENVIRONMENT_NAME
conda install -n ENVIRONMENT_NAME -c conda-forge rdkit pandas numpy biopython jinja2 griddataformats parmed
conda install -n ENVIRONMENT_NAME -c conda-forge openbabel==2.4.1
```


### Modification of `setting/initialize.sh`

Please modify `setting/initialize.sh` to correctly use all executables and GPUs.

### Test execution
If all settings goes well, below commands run MD calculations.
```
cd /PATH/TO/exprorer/cmd_calculation              # MUST RUN THEM AT THIS DIRECTORY
bash msmd_preparation.sh example/msmd_config.sh   # system prepration
bash msmd_mdrun.sh example/msmd_config.sh         # MD execution
```
The commands will execute 20 runs of CMD with Gromacs, each of which includes 40 ns of production run as a default.

### construct spatial probability map (PMAP)
```
python script/combine_pmaps.py -i [dx files of all runs] -o pmap.dx \
       --mode probability -m [reference protein structure].pdb
```
The command will construct PMAP in accordance with input map files.

## How to modify config files

### To change the number of runs / the length of each production run

1. Modify `msmd_config.sh`

### To do CMD with a different protein

1. Prepare the protein `.pdb` file
  - Residue IDs must be 1-origin and must not any jump of residue IDs.
  - The protein should be centered.
2. Make new protein config file (example: `example/protein.conf`). You have to write down an absolute path to protein PDB file and the residue numbers which make disulfide bonds. 
3．Make new MSMD config file like `example/msmd_config.sh`. Please make sure to use new protein config file.

### To do CMD with a different probe (cosolvent)

1. Prepare probe `.mol2` and `.pdb` files
2. Make new probe config file (example: `example/probe.conf`). You have to write down absolute pathes to `.mol2` and `.pdb` files as well as residue name of the probe.
3．Make new MSMD config file like `example/msmd_config.sh`. Please make sure to use new probe config file.

## Calculation cost

- 3.5 hrs for 40 ns simulation (with 37,763 atoms of example protocol) on Cygnus, Tsukuba University (Xeon Gold 6126 + Tesla V100)
  - 70 GPU hrs for 20 runs of 40 ns simulation