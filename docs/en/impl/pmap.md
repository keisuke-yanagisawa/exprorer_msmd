# PMAP (Probability Distribution Map)

PMAP (spatial Probability distribution MAP) is a three-dimensional probability density map representing the spatial distribution of probe molecules.

## Table of Contents
1. [PMAP Basics](#pmap-basics)
2. [Calculation Method](#calculation-method)
   - [Basic Procedure](#1-basic-procedure)
   - [Grid Setup](#2-grid-setup)
   - [Atom Selection](#3-atom-selection)
   - [Value Conversion](#4-value-conversion)
3. [Implementation Details](#implementation-details)
   - [Single Simulation Processing](#1-single-simulation-processing-gen_pmap-function)
   - [Multiple Simulation Integration](#2-multiple-simulation-integration)

Related Documents:
- [MSMD Implementation](msmd.md)
- [Basic Usage](../user_guide/basic.md)
- [Advanced Usage](../user_guide/advanced.md)

## PMAP Basics

PMAP represents the spatial probability distribution of probe molecules observed in mixed-solvent molecular dynamics (MSMD) simulations.

This information can be utilized both qualitatively and quantitatively:

- Regions where probes are likely to exist are likely binding sites for compounds (qualitative, hotspot analysis)
- Regions with high probe existence probability likely indicate stronger binding (quantitative, binding affinity prediction)

## Calculation Method

### 1. Basic Procedure

1. Division of simulation box into 3D grid
2. Counting probe atom frequencies at each grid point
3. Converting frequencies to probabilities

### 2. Grid Setup

Grid creation is performed using cpptraj from Amber Tools.

- Grid center: Protein center of mass (set as origin)
- Default box size: 80 Å × 80 Å × 80 Å
- Grid spacing: 1 Å

This creates an 80 × 80 × 80 grid with 1 Å spacing,
counting the **frequency** of probe atoms at each grid point.

### 3. Atom Selection

Atoms used for grid creation are selected by selectors specified in the yaml file.
Here are three examples:

1. Heavy atoms only
   - Selector: `(!@VIS)&(!@H*)`
   - Excludes hydrogen atoms and virtual atoms
2. All atoms
   - Selector: `(!@VIS)`
   - Excludes only virtual atoms
3. Probe center of mass
   - Selector: `@VIS`
   - Selects only probe virtual atoms
     - Since probe virtual atoms are placed at the center of mass of each probe, selecting only virtual atoms effectively selects probe centers of mass.

### 4. Value Conversion

As mentioned above, the grid created so far counts the **frequency** of probe atoms.
Naturally, these counts vary with simulation length (more precisely, the number of snapshots),
so they need to be converted to **probabilities**. Two types of normalization are provided.
The default is **total normalization** (although snapshot normalization may often be more appropriate).

1. Total normalization
   $P(r) = N(r) / \sum N(r)$

   - $N(r)$: Frequency at position $r \in \mathbb{R}^3$
   - Sum of probabilities across all grid points equals 1

2. Snapshot normalization
   $P(r) = N(r) / N_{\mathrm{frames}}$

   - $N_{\mathrm{frames}}$: Number of snapshots used
   - Directly represents existence probability at each position

Additionally, a method is provided to estimate the interaction energy between protein and probe molecules (grid free energy; GFE) using these probabilities.

$$\mathrm{GFE}(r) = -RT \ln(P(r)/P_{\mathrm{bulk}})$$

- $R$: Gas constant
- $T$: Temperature
- $P_{\mathrm{bulk}}$: Bulk probability (average probability across the entire space)

## Implementation Details

### 1. Single Simulation Processing (`gen_pmap()` function)

#### Frequency Grid Calculation with cpptraj

Grid calculation using cpptraj follows these steps.
The template for this execution can be found in `script/utilities/executables/template/cpptraj_pmap.in`.

1. Trajectory preprocessing
   Loads the trajectory and performs alignment with the initial structure using backbone carbon atoms.
   ```
   trajin [project_name].xtc [start_frame] [end_frame]
   reference [project_name].pdb
   rms reference @CA
   ```

2. Frequency grid creation
   Creates a grid based on the specified atom selector using cpptraj's grid command.

   ```
   grid [output_map_file_path].dx
     [X_size(Å)] [X_spacing(Å)] [Y_size(Å)] [Y_spacing(Å)] [Z_size(Å)] [Z_spacing(Å)]
     gridcenter 0.0 0.0 0.0 [atom_selector]
   ```
   - Default size: 80×80×80 Å
   - Default spacing: 1×1×1 Å

#### PMAP File Generation

Reads the frequency grid created by cpptraj and generates PMAP.

1. Conversion to probabilities
   - For total normalization:
     ```python
     total_count = np.sum(frequency_map)
     probability_map = frequency_map / total_count
     ```
   - For snapshot normalization:
     ```python
     probability_map = frequency_map / n_frames
     ```

2. GFE calculation (optional) `convert_to_gfe()`
   ```python
   R = 0.001987  # kcal/mol/K
   T = 300       # K
   bulk_prob = np.mean(probability_map)
   gfe_map = -R * T * np.log(probability_map / bulk_prob)
   ```

   Note that `gfe_map` may contain infinite values (this occurs at coordinates where probability_map is 0),
   so if any value in `gfe_map` is greater than 3.0 kcal/mol, it is clipped to 3.0 kcal/mol.

### 2. Multiple Simulation Integration

under construction

Note: All maps are output in OpenDX format.