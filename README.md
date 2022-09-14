# EXPRORER_MSMD

MSMD (mixed-solvent molecular dynamics) engine and analysis tools

## Requirements

Below requirements are written as `.devcontainer/Dockerfile`. You may use docker container.

- python 3 (`$PYTHON`)
  - Several modules are needed (an example of installation is shown below).
- AmberTools 20 (`$TLEAP`, `$CPPTRAJ`)
  - It is needed for use of `tleap` and `cpptraj`
- Gromacs 2021.5 (`$GMX`)
- packmol 18.169 (`$PACKMOL`)

## Preparation

### Preparation of container environment

- Docker: The container can be built with `Dockerfile`.
- Singularity: Firstly, convert `Dockerfile` to a Singularity Recipe with Singularity Python. Then, build the container with Singularity.
It can be used with Docker. Below instruction is for user who will not use Docker.

### Test execution
If all settings goes well, below commands run MD calculations.
```
cd /PATH/TO/exprorer_msmd
./exprorer_msmd   example/example_protocol.yaml # prepare a system & execute MSMD simulation
./protein_hotspot example/example_protocol.yaml # construct a probability map on protein surface
./probe_profile   example/example_protocol.yaml # construct an interaction profile of a probe
```

## How to modify config files

### To change the number of runs / the length of each production run

1. Modify the protocol yaml file.

### To do CMD with a different protein

1. Prepare the protein `.pdb` file
  - Residue IDs must be 1-origin and must not any jump of residue IDs.
2．Modify or make a protocol file with settings of `input` - `protein`. 

### To do CMD with a different probe (cosolvent)

1. Prepare probe `.mol2` and `.pdb` files
2．Modify or make a protocol file with settings of `input` - `probe`.

## Calculation cost

- 3.5 hrs for 40 ns simulation (with 37,763 atoms of example protocol) on Cygnus, Tsukuba University (Xeon Gold 6126 + Tesla V100)
  - 70 GPU hrs for 20 runs of 40 ns simulation