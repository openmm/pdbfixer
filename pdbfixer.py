"""
pdbfixer.py: Fixes problems in PDB files

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2013 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
__author__ = "Peter Eastman"
__version__ = "1.0"

import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
from simtk.openmm.app.element import hydrogen, oxygen
import numpy as np
import numpy.linalg as lin
import sys
import os
import os.path
import math

def overlayPoints(points1, points2):
    """Given two sets of points, determine the translation and rotation that matches them as closely as possible.
    
    This is based on W. Kabsch, Acta Cryst., A34, pp. 828-829 (1978)."""
    
    if len(points1) == 0:
        return (Vec3(0, 0, 0), np.identity(3), Vec3(0, 0, 0))
    if len(points1) == 1:
        return (points1[0], np.identity(3), -points2[0])
    
    # Compute centroids.
    
    center1 = unit.sum(points1)/float(len(points1))
    center2 = unit.sum(points2)/float(len(points2))
    
    # Compute R matrix.
    
    R = np.zeros((3, 3))
    for p1, p2 in zip(points1, points2):
        x = p1-center1
        y = p2-center2
        for i in range(3):
            for j in range(3):
                R[i][j] += y[i]*x[j]
    
    # Use an SVD to compute the rotation matrix.
    
    (u, s, v) = lin.svd(R)
    return (-1*center2, np.dot(u, v).transpose(), center1)

# Load the PDB file.

pdb = app.PDBFile(sys.argv[1])

# Load the templates.

templates = {}
for file in os.listdir('templates'):
    templatePdb = app.PDBFile(os.path.join('templates', file))
    name = templatePdb.topology.residues().next().name
    templates[name] = templatePdb

# Loop over residues to see which ones have missing heavy atoms.

missingAtoms = {}
for residue in pdb.topology.residues():
    if residue.name in templates:
        template = templates[residue.name]
        atomNames = set(atom.name for atom in residue.atoms())
        missing = []
        for atom in template.topology.atoms():
            if atom.name not in atomNames:
                missing.append(atom)
        if len(missing) > 0:
            missingAtoms[residue] = missing

# Create the new Topology.

newTopology = app.Topology()
newPositions = []*unit.nanometer
existingAtomMap = {}
addedAtomMap = {}
addedOXT = []
for chain in pdb.topology.chains():
    newChain = newTopology.addChain()
    chainResidues = list(chain.residues())
    for residue in chain.residues():
        newResidue = newTopology.addResidue(residue.name, newChain)
        
        # Add the existing heavy atoms.
        
        for atom in residue.atoms():
            if atom.element is not None and atom.element != hydrogen:
                newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                existingAtomMap[atom] = newAtom
                newPositions.append(pdb.positions[atom.index])
        if residue in missingAtoms:
            
            # Find corresponding atoms in the residue and the template.
            
            template = templates[residue.name]
            atomPositions = dict((atom.name, pdb.positions[atom.index]) for atom in residue.atoms())
            points1 = []
            points2 = []
            for atom in template.topology.atoms():
                if atom.name in atomPositions:
                    points1.append(atomPositions[atom.name].value_in_unit(unit.nanometer))
                    points2.append(template.positions[atom.index].value_in_unit(unit.nanometer))
            
            # Compute the optimal transform to overlay them.
            
            (translate2, rotate, translate1) = overlayPoints(points1, points2)
            
            # Add the missing atoms.
            
            addedAtomMap[residue] = {}
            for atom in missingAtoms[residue]:
                newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                addedAtomMap[residue][atom] = newAtom
                templatePosition = template.positions[atom.index].value_in_unit(unit.nanometer)
                newPositions.append((mm.Vec3(*np.dot(rotate, templatePosition+translate2))+translate1)*unit.nanometer)
        
        # If a terminal OXT is missing, add it.
        
        if residue == chainResidues[-1] and residue.name in templates:
            atomPositions = dict((atom.name, pdb.positions[atom.index].value_in_unit(unit.nanometer)) for atom in residue.atoms())
            if 'OXT' not in atomPositions and all(name in atomPositions for name in ['C', 'O', 'CA']):
                newAtom = newTopology.addAtom('OXT', oxygen, newResidue)
                addedOXT.append(newAtom)
                d_ca_o = atomPositions['O']-atomPositions['CA']
                d_ca_c = atomPositions['C']-atomPositions['CA']
                d_ca_c /= unit.sqrt(unit.dot(d_ca_c, d_ca_c))
                v = d_ca_o - d_ca_c*unit.dot(d_ca_c, d_ca_o)
                newPositions.append((atomPositions['O']+2*v)*unit.nanometer)

# Add bonds from the original Topology.

for atom1, atom2 in pdb.topology.bonds():
    if atom1 in existingAtomMap and atom2 in existingAtomMap:
        newTopology.addBond(existingAtomMap[atom1], existingAtomMap[atom2])

# Add bonds that connect to new atoms.

for residue in missingAtoms:
    template = templates[residue.name]
    atomsByName = dict((atom.name, atom) for atom in residue.atoms())
    addedAtoms = addedAtomMap[residue]
    for atom1, atom2 in template.topology.bonds():
        if atom1 in addedAtoms or atom2 in addedAtoms:
            if atom1 in addedAtoms:
                bondAtom1 = addedAtoms[atom1]
            else:
                bondAtom1 = existingAtomMap[atomsByName[atom1.name]]
            if atom2 in addedAtoms:
                bondAtom2 = addedAtoms[atom2]
            else:
                bondAtom2 = existingAtomMap[atomsByName[atom2.name]]
            newTopology.addBond(bondAtom1, bondAtom2)
for atom1 in addedOXT:
    atom2 = [atom for atom in atom1.residue.atoms() if atom.name == 'C'][0]
    newTopology.addBond(atom1, atom2)

app.PDBFile.writeFile(newTopology, newPositions, open('output.pdb', 'w'))

forcefield = app.ForceField('soft.xml')
forcefield.createSystem(newTopology)