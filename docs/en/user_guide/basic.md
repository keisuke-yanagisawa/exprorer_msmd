# Basic Usage

## Introduction

EXPRORER_MSMD is a tool for running and analyzing Mixed-Solvent Molecular Dynamics (MSMD) simulations.

## Environment Setup

This repository includes Docker files that allow you to easily create an execution environment with the following tools:

- Python 3
- AmberTools 20
- Gromacs 2021.5
- Packmol 18.169

## Basic Usage Steps

### 1. Protocol File Preparation

Simulation settings are defined in a YAML file:

```yaml
general:
  name: TEST_PROJECT         # Project name
  workdir: ./PATH/TO/WORKDIR # Output directory

input:
  protein:
    pdb: protein.pdb       # Protein structure
  probe:
    cid: PROBE            # Probe molecule name
    molar: 0.25          # Concentration (mol/L)
```

Please also refer to `example/example_protocol.yaml` for YAML file descriptions.

### 2. Running the Simulation

```bash
./exprorer_msmd protocol.yaml
```

This will automatically perform system construction, simulation, and PMAP creation.

### 3. Analyzing Results

#### Checking Trajectories

The simulation results are automatically converted to PDB files and
saved as `./PATH/TO/WORKDIR/system*/[project_name]_woWAT_10ps.pdb`.
Note that `*` represents the simulation ID at runtime, and when there are multiple independent trials, each will have its own ID.

These trajectories can be opened with PyMOL to observe the movements of proteins and probes during the simulation.

#### Hotspot Analysis
```bash
./protein_hotspot protocol.yaml
```

Under construction

#### Probe Environment Analysis
```bash
./probe_profile protocol.yaml
```

Under construction