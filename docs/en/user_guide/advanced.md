# Advanced Usage

## Table of Contents
1. [Introduction](#introduction)
2. [Customizing Simulation Settings](#customizing-simulation-settings)
   - [Controlling Simulation Steps](#controlling-simulation-steps)
   - [Independent Trial Settings](#independent-trial-settings)
   - [Adjusting Simulation Parameters](#adjusting-simulation-parameters)
3. [Customizing Analysis Settings](#customizing-analysis-settings)
   - [PMAP Calculation Settings](#pmap-calculation-settings)
   - [Inverse MSMD Related Settings](#inverse-msmd-related-settings)

## Introduction

This document explains the detailed configuration options for EXPRORER_MSMD.
Please refer to this after understanding [Basic Usage](basic.md) when you want to have finer control over simulations and analyses.

## Customizing Simulation Settings

We provide options to control simulation behavior in detail.

### Controlling Simulation Steps

Options are available to execute each processing phase individually.
This is useful when you want to reanalyze with different parameters or try analysis with different conditions.

```bash
# Skip preprocessing
# - When using an existing system
# - When system construction is complete but simulation needs to be rerun
./exprorer_msmd protocol.yaml --skip-preprocess

# Skip simulation
# - When performing different analyses on existing trajectories
# - When reanalyzing with modified analysis parameters
./exprorer_msmd protocol.yaml --skip-simulation

# Skip postprocessing
# - When only running simulations
# - When analysis will be performed separately
./exprorer_msmd protocol.yaml --skip-postprocess
```

### Independent Trial Settings

You can run multiple independent simulations to increase reliability.

```yaml
general:
  iter_index: 0,1,2  # Execute 3 independent trials
  # Notation:
  # "1-3" => 1,2,3
  # "5-9:2" => 5,7,9
  # "1-3,5-9:2" => 1,2,3,5,7,9
```

### Adjusting Simulation Parameters

You can set parameters to control the physical conditions of the simulation.

```yaml
exprorer_msmd:
  general:
    dt: 0.002          # Time step (ps)
    temperature: 300   # Temperature (K)
    pressure: 1.0      # Pressure (bar)
    pbc: xyz          # Periodic boundary conditions

  sequence:
    - name: pr        # Production run
      type: production
      nsteps: 20000000  # Number of steps (40 ns in this example)
      nstxtcout: 5000   # Output frequency (every 10 ps in this example)
```

## Customizing Analysis Settings

We provide options to control how simulation results are analyzed.

### PMAP Calculation Settings

You can control the calculation method of Probability MAP (PMAP) in detail.

```yaml
map:
  type: pmap
  snapshot: 2001-4001:1  # Snapshot range to use: use data after equilibration
  maps:
    - suffix: nVH        # Map using heavy atoms only
      selector: (!@VIS)&(!@H*)
    - suffix: nV         # Map using all atoms
      selector: (!@VIS)
  map_size: 80           # Map size (Å): specify a size large enough to contain the entire system
  normalization: total   # Normalization method: total, snapshot, or GFE can be specified
```

### Inverse MSMD Related Settings

You can configure settings for probe molecule environment analysis.

```yaml
probe_profile:
  resenv:
    map: nVH            # Map to use: typically use heavy-atoms-only map
    threshold: 0.001    # Probability threshold: lower values detect wider range of interactions
  profile:
    types:
      - name: anion     # Interaction with anionic residues
        atoms:          # Specify atoms for interaction evaluation
          - ["ASP", " CB "]
          - ["GLU", " CB "]
```

For detailed explanations and recommended values for each setting, please also refer to `example/example_protocol.yaml`.