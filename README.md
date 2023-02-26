[![PDBFixer Continuous Integration Workflow](https://github.com/openmm/pdbfixer/actions/workflows/CI.yml/badge.svg)](https://github.com/openmm/pdbfixer/actions/workflows/CI.yml)

PDBFixer
========

PDBFixer is an easy to use application for fixing problems in Protein Data Bank files in preparation for simulating them.  It can automatically fix the following problems:

- Add missing heavy atoms.
- Add missing hydrogen atoms.
- Build missing loops.
- Convert non-standard residues to their standard equivalents.
- Select a single position for atoms with multiple alternate positions listed.
- Delete unwanted chains from the model.
- Delete unwanted heterogens.
- Build a water box for explicit solvent simulations.

See our [manual](https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html)

## Installation

PDBFixer can be installed with conda or mamba.

```
conda install -c conda-forge pdbfixer
```

Alternatively you can install from source, as described in the manual.