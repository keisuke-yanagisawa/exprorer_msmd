last modified: 19-May-2023

# EXPRORER_MSMD

Repository that combines the Mixed-Solvent Molecular Dynamics (MSMD) simulation engine and analysis tools.

## System Description

A system for performing MSMD using GROMACS automatically.

## Usage Guide

### Environment Setup

This GitHub repository includes Docker files, which can be used to easily create an executable environment.

The notable applications used inside are:

- Python 3 (`$PYTHON`)
  - Refer to `.devcontainer/Dockerfile` for the libraries used.
- AmberTools 20 (`$TLEAP`, `$CPPTRAJ`)
- Gromacs 2021.5 (`$GMX`)
- Packmol 18.169 (`$PACKMOL`)

### Performing MSMD Simulations (`exprorer_msmd`)

Prepare a YAML file that defines the protein, mixed-solvent (probe) molecules, and simulation protocol (e.g., `example/example_protocol.yaml`). By executing the following command, the system will automatically perform the system setup, simulation, and create a spatial probability distribution map (PMAP).

```
./exprorer_msmd example/example_protocol.yaml
```


It should be noted that the simulation procedure consists of three parts:

- Preprocessing to create the MSMD system (**preprocess**)
- Simulation using GROMACS (**simulation**)
- Postprocessing to generate PMAP files from simulation results (**postprocess**)

Each part can be skipped using the flags `--skip-preprocess`, `--skip-simulation`, and `--skip-postprocess`, respectively.

### Analysis of MSMD Simulation Results

#### Checking the Simulation Trajectory

When executing `./exprorer_msmd`, `system` folders will be created in the output directory for each independent run. Within this folder, several PDB files such as `[project_name]_woWAT_10ps.pdb` will be generated, which represent the MSMD trajectory. You can use PyMOL to open these files and visualize the movements of the protein and probe molecules during the simulation (water molecules are excluded).

#### Protein Hotspot Exploration (`protein_hotspot`)

This step involves exploring the regions on the protein surface where the probe molecules are most likely to be found, known as protein hotspots.

```
./protein_hotspot example/example_protocol.yaml
```

This command will generate voxel files in OpenDX format, such as `maxPMAP_[project_name]_nV.dx`, in the root directory of the output. These voxel files are created to align with the input protein structure. You can load them in PyMOL and visualize them using the `isomesh` command.

The resulting visualization may resemble the following image (showing multiple probe calculation results overlaid):

![Protein Hotspot Exploration](https://i.imgur.com/bzxz0K6.png)

EXPRORER[^1] uses these results to calculate the similarity between voxels and between probes.

[^1]: **Keisuke Yanagisawa**, Yoshitaka Moriwaki, Tohru Terada, Kentaro Shimizu. "EXPRORER: Rational Cosolvent Set Construction Method for Cosolvent Molecular Dynamics Using Large-Scale Computation", *Journal of Chemical Information and Modeling*, **61**: 2744-2753, 2021/06. DOI: [10.1021/acs.jcim.1c00134](https://doi.org/10.1021/acs.jcim.1c00134)

#### Obtaining the Residue Environment around Probe Molecules[^2] (`probe_profile`)

In contrast to the hotspot analysis mentioned earlier, this step involves visualizing the residues that are most likely to be found around the probe molecules.

```
./protein_hotspot example/example_protocol.yaml
```

This command will generate eight voxel files in OpenDX format, such as `[project_name]_[probe_name]_mesh_anion.dx`, in the root directory of the output. These voxel files are created to align with the `alignedresenv_[project_name].pdb` file, and you can visualize them in PyMOL using the `isomesh` command.

The resulting visualization may resemble the following image:
![Residue Environment Exploration](https://i.imgur.com/4QIZxhW.png)

[^2]: **Keisuke Yanagisawa**, Ryunosuke Yoshino, Genki Kudo, Takatsugu Hirokawa. "Inverse Mixed-Solvent Molecular Dynamics for Visualization of the Residue Interaction Profile of Molecular Probes", *International Journal of Molecular Sciences*, **23**: 4749, 2022/04. DOI: [10.3390/ijms23094749](https://doi.org/10.3390/ijms23094749)

## Editing the YAML File

All configurations are specified in the YAML file.
(Under construction)
