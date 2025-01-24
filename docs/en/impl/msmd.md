# Implementation Details

Mixed-Solvent Molecular Dynamics (MSMD) is a molecular dynamics simulation method for analyzing interactions between proteins and multiple solvent molecules (probes).

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