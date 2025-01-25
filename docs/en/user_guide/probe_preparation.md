# Probe Molecule Preparation

Probe molecules are prepared following these steps:

1. Initial 3D structure generation
2. Probe molecule structure optimization
3. Electron density calculation
4. Partial charge assignment using RESP method

Note that probe molecules should ideally **have no internal degrees of freedom**. Especially when calculating the probability of probe existence on protein surfaces, the conformational state of the probe molecule itself has a significant impact. By eliminating internal degrees of freedom in probe molecules, we can prevent reduced sampling efficiency related to probe molecule conformations.

If you need to use probe molecules with internal degrees of freedom, make sure to increase the total amount of MSMD simulation and verify that the results have sufficiently converged.

As an example, `./example/A11.mol2` and `./example/A11.pdb` are registered as isopropanol probes with no internal degrees of freedom.

## 1. Initial 3D Structure Generation

Use tools like RDKit or OpenBabel to generate the initial 3D structure of probe molecules.
Since quantum chemical calculations will be performed starting from this initial structure, carefully verify that the generated initial 3D structure meets your expectations (e.g., for probe molecules with chiral centers or internal degrees of freedom).

```sh
# Example of generating initial 3D structure for isopropanol using OpenBabel
# (Prepare SMILES notation for isopropanol in advance)
obabel -ismi isopropanol.smi -omol2 -Oisopropanol.mol2 --gen3D
```

## 2. Probe Molecule Structure Optimization

Next, perform structure optimization of the probe molecule using quantum chemical calculations. Insufficient optimization at this stage may affect the partial charge assignment in the next step.

While there are various quantum chemistry software packages like Gaussian and GAMESS, here's an example using Gaussian 16.

First, prepare the Gaussian 16 input file. exprorer_msmd uses the B3LYP/6-31G(d) level of theory.

`./example/A11_opt.gjf`
```
%chk=opt.chk
%mem=10GB
# p opt=(maxcycle=999, maxstep=5, tight) b3lyp/6-31g(d) scf=(qc,tight)
opt

0 1
0  1
C           1.06720        -0.04599        -0.01121 # Atomic coordinates follow
C           2.58785        -0.07121        -0.02775
...
```

Input this file to Gaussian 16:

```sh
g16 < A11_opt.gjf > A11_opt.log
```

The execution log will be output to `A11_opt.log`. If you find the string `Normal termination` in it, the calculation completed successfully.

## 3. Electron Density Calculation

Next, calculate the electron density for the optimized structure. Use Gaussian 16's `pop=mk` option for this calculation. Note that this calculation should be performed at the HF/6-31g(d) level, **not** at the B3LYP/6-31g(d) level.

`./example/A11_elec.gjf`
```
%chk=hf.chk
%mem=10GB
# p HF/6-31g(d) iop(6/33=2,6/41=10,6/42=17) pop=(mk,readradii) scf=tight
hf

0 1
0  1
C           1.06720        -0.04599        -0.01121 # Atomic coordinates follow
C           2.58785        -0.07121        -0.02775 # Use post-optimization coordinates (note: these coordinates are pre-optimization)
...
```

```sh
g16 < A11_elec.gjf > A11_elec.log
```

This will calculate the electron density and record the results in the log file.

## 4. Partial Charge Assignment using RESP Method

Finally, assign partial charges using the RESP method. This is done using `antechamber`.

```sh
# Run antechamber twice to output both mol2 and pdb files
# -nc specifies the total molecular charge. Isopropanol is neutral, so 0
# -rn specifies the residue name. Choose a name that won't conflict with other molecules
antechamber \
  -fi gout -i A11_elec.log \
  -fo mol2 -o A11_resp.mol2 \
  -c resp \
  -nc 0 \
  -at gaff2 \
  -pf y \
  -rn A11

antechamber \
  -fi gout -i A11_elec.log \
  -fo pdb -o A11_resp.pdb \
  -c resp \
  -nc 0 \
  -at gaff2 \
  -pf y \
  -rn A11
```

Note that when outputting mol2 files with `antechamber`, sometimes the sum of partial charges may not equal zero. This will cause errors during MSMD simulation execution, so in such cases, manually adjust the charges to ensure their sum equals zero.