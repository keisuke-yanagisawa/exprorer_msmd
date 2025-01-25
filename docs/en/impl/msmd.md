# Implementation Details

Mixed-Solvent Molecular Dynamics (MSMD) is a molecular dynamics simulation method for analyzing interactions between proteins and multiple solvent molecules (probes).

## Table of Contents
1. [MSMD Basics](#msmd-basics)
   - [Probe Molecules](#probe-molecules)
   - [Virtual Repulsive Forces Between Probe Molecules](#virtual-repulsive-forces-between-probe-molecules)
2. [MSMD Simulation Protocol](#msmd-simulation-protocol)
   - [Energy Minimization](#energy-minimization)
   - [System Heating](#system-heating)
   - [Gradual Release of Position Restraints](#gradual-release-of-position-restraints)
   - [Production Run](#production-run)
3. [Implementation Details](#implementation-details-1)
   - [Preprocessing](#1-preprocessing-preprocess)
   - [Simulation Execution](#2-simulation-execution-execute_single_simulation)
   - [Postprocessing](#3-postprocessing-postprocess)

Related Documents:
- [PMAP Details](pmap.md)
- [Basic Usage](../user_guide/basic.md)
- [Probe Preparation](../user_guide/probe_preparation.md)

## MSMD Basics

MSMD simulations are characterized by:

1. Multiple probe molecules present in the system
2. Virtual repulsive forces between probe molecules
3. Multiple short simulations with different initial structures

### Probe Molecules

Probe molecules are analytical molecules in MSMD, typically very small organic compounds such as benzene or isopropanol. These probe molecules interact with proteins, enabling various analyses such as promoting protein structural changes (cryptic binding site exploration) and analyzing molecular binding on protein surfaces (hotspot analysis, binding affinity prediction).

### Virtual Repulsive Forces Between Probe Molecules

When using hydrophobic probe molecules, they tend to aggregate. To solve this problem, virtual repulsive forces are introduced between probe molecules.

The virtual repulsive force between probe molecules is expressed using the Lennard-Jones potential:

1. Place a virtual atom VIS at the center of mass of each probe molecule
2. Apply Lennard-Jones potential between virtual atoms VIS

$U_{LJ}(r) = 4\epsilon_{LJ} \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]$

- $\epsilon_{LJ} = 10^{-6} \mathrm{kcal/mol}$
- $R_{min} = 20.0 \mathrm{\AA}$ ($\sigma = R_{min}/2^{1/6}$)

These parameters were determined by manually searching for values that well reproduce the carbon-carbon distances between probe molecules shown in Guvench et al.'s research[^1].

[^1]: Olgun Guvench and Alexander D. MacKerell Jr. "Computational Fragment-Based Binding Site Identification by Ligand Competitive Saturation", *PLoS Computational Biology*, **5**: e1000435, 2009. DOI: [10.1371/journal.pcbi.1000435](https://doi.org/10.1371/journal.pcbi.1000435)

## MSMD Simulation Protocol

In MSMD, the following simulation steps are performed multiple times with varying initial system structures.

### Energy Minimization
To partially resolve strain in the system, steepest descent energy minimization is performed.
First, minimization is performed with position restraints, followed by minimization without restraints.

### System Heating
The system is heated to the target temperature.

- NVT ensemble
- Position restraints on protein
- Gradual temperature increase

### Gradual Release of Position Restraints
From this point, as an equilibration step, position restraints are gradually released.

1. Start with strong restraints (1000 kcal/mol/Å²)
2. Gradually weaken restraints
3. Finally remove all restraints

The restraint strength changes in the following stages:

1. 1000 kcal/mol/Å²
2. 500 kcal/mol/Å²
3. 200 kcal/mol/Å²
4. 100 kcal/mol/Å²
5. 50 kcal/mol/Å²
6. 20 kcal/mol/Å²
7. 10 kcal/mol/Å²
8. 0 kcal/mol/Å²

These are performed under constant pressure (NPT ensemble).

### Production Run
The final simulation is executed.

- NPT ensemble
- No position restraints
- Fixed simulation time (standard: 40 ns)

## Implementation Details

`./exprorer_msmd` consists of three main processing stages:

1. Preprocessing `preprocess()`
2. Simulation Execution `execute_single_simulation()`
3. Postprocessing `postprocess()`

For the following explanation, we'll use an example where the simulation iteration ID is `42`.

### 1. Preprocessing `preprocess()`

Preprocessing involves building the system necessary for MSMD simulation. The final output data is converted to a format that can be handled by GROMACS.
The execution results are saved in the `PATH/TO/OUTPUT/DIR/system42/prep` directory.

#### 1-1. Probe Preparation

- Reference the probe's mol2 file and use the `parmchk2` command to generate an frcmod file that supplements missing force field parameters.

#### 1-2. Protein Preparation

- Read the input PDB file and remove ANISOU records (anisotropic temperature factor information) and C-terminal oxygen atoms (OXT).
- Determine the system size from the protein structure.
   - Create a system containing only protein and water molecules using `tleap`, and record the length of the longest edge of this system.
   - Use a cube based on this length as the size of the MSMD simulation system.
     - A cubic region was chosen to prevent interference between proteins across boundaries, especially in cases where the protein has an elongated shape.

#### 1-3. MSMD System Construction

- Create a system containing only protein and probe molecules.
  - Calculate the number of probe molecules that can be placed in the system based on the specified probe molar concentration and the system size determined in 1-2.
  - Use `packmol` to place the protein and probe molecules within the system. The initial placement of probe molecules is randomized for each simulation trial using the simulation ID as a seed value.
- Add water molecules and ions
  - Use `tleap` to fill the system with water molecules and add ions to neutralize the system charge.
    - Na+ and Cl- ions are added only as minimally necessary, which may result in cases where neither Na+ nor Cl- is added.
- Convert topology files
  - Convert Amber format topology files (.parm7/.rst7) to GROMACS format (.top/.gro).
- Add virtual atoms
   - Place virtual atoms VIS at the center of mass of probe molecules.
   - Add VIS coordinates to the gro file
   - Set up `[virtual_sitesn]` in the top file to fix VIS positions.
   - Add non-bonded interactions between VIS-VIS pairs in the top file to apply Lennard-Jones potential.
- Set up position restraints
   - Add entries to the top file to enable 8 stages of position restraints for protein heavy atoms.

### 2. Simulation Execution `execute_single_simulation()`

Using the system created in preprocessing, MSMD simulation is executed. GROMACS is used as the simulation engine.
The execution results are saved in the `PATH/TO/OUTPUT/DIR/system42/simulation` directory.

- Generate index file
   - Define atom groups using `gmx make_ndx`. exprorer_msmd uses two groups: `Water` and `Non-Water`.
- Create mdp files
   - Automatically generate mdp files for each step (minimization, heating, equilibration, production) based on YAML settings.
      - Templates for each step exist in `./script/template/production.mdp` etc.
- Execute simulation
   - Automatically create and execute `mdrun.sh` that runs each simulation step sequentially.
   - A symbolic link named `[project_name].xtc` is created for the trajectory file obtained in the final step (typically pr (production run)).

### 3. Postprocessing `postprocess()`

For details about this step, please refer to [PMAP Details](pmap.md).