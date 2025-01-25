# Basic Usage

## Table of Contents
1. [Introduction](#introduction)
2. [Environment Setup](#environment-setup)
3. [Basic Usage Steps](#basic-usage-steps)
   - [Protocol File Preparation](#1-protocol-file-preparation)
   - [Running the Simulation](#2-running-the-simulation)
   - [Analyzing Results](#3-analyzing-results)

## Introduction

EXPRORER_MSMD is a tool for running and analyzing Mixed-Solvent Molecular Dynamics (MSMD) simulations.
This tool enables the following types of analysis:

- Distribution analysis of co-solvent molecules on protein surfaces
- Protein hotspot identification
- Probe molecule environment analysis

## Environment Setup

This repository includes Docker files that allow you to easily create an execution environment with the following tools:

- Python 3
- AmberTools 20
- Gromacs 2021.5
- Packmol 18.169

## Basic Usage Steps

### 1. Protocol File Preparation

Simulation settings are defined in a YAML file. Below is a partial example of the settings.
For a complete example, please refer to `example/example_protocol.yaml`.
For more detailed configuration options, please refer to [Advanced Usage](advanced.md).

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

### 2. Running the Simulation

```bash
./exprorer_msmd protocol.yaml
```

This will automatically perform the following processes:

1. System construction
2. MD simulation execution
3. PMAP creation

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

## Next Steps

- For more detailed configuration methods, please refer to [Advanced Usage](advanced.md).
- For information on preparing probe molecules, please refer to [Probe Preparation](probe_preparation.md).