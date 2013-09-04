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

import ui
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

substitutions = {
    '2AS':'ASP', '3AH':'HIS', '5HP':'GLU', 'ACL':'ARG', 'AGM':'ARG', 'AIB':'ALA', 'ALM':'ALA', 'ALO':'THR', 'ALY':'LYS', 'ARM':'ARG',
    'ASA':'ASP', 'ASB':'ASP', 'ASK':'ASP', 'ASL':'ASP', 'ASQ':'ASP', 'AYA':'ALA', 'BCS':'CYS', 'BHD':'ASP', 'BMT':'THR', 'BNN':'ALA',
    'BUC':'CYS', 'BUG':'LEU', 'C5C':'CYS', 'C6C':'CYS', 'CCS':'CYS', 'CEA':'CYS', 'CGU':'GLU', 'CHG':'ALA', 'CLE':'LEU', 'CME':'CYS',
    'CSD':'ALA', 'CSO':'CYS', 'CSP':'CYS', 'CSS':'CYS', 'CSW':'CYS', 'CSX':'CYS', 'CXM':'MET', 'CY1':'CYS', 'CY3':'CYS', 'CYG':'CYS',
    'CYM':'CYS', 'CYQ':'CYS', 'DAH':'PHE', 'DAL':'ALA', 'DAR':'ARG', 'DAS':'ASP', 'DCY':'CYS', 'DGL':'GLU', 'DGN':'GLN', 'DHA':'ALA',
    'DHI':'HIS', 'DIL':'ILE', 'DIV':'VAL', 'DLE':'LEU', 'DLY':'LYS', 'DNP':'ALA', 'DPN':'PHE', 'DPR':'PRO', 'DSN':'SER', 'DSP':'ASP',
    'DTH':'THR', 'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL', 'EFC':'CYS', 'FLA':'ALA', 'FME':'MET', 'GGL':'GLU', 'GL3':'GLY', 'GLZ':'GLY',
    'GMA':'GLU', 'GSC':'GLY', 'HAC':'ALA', 'HAR':'ARG', 'HIC':'HIS', 'HIP':'HIS', 'HMR':'ARG', 'HPQ':'PHE', 'HTR':'TRP', 'HYP':'PRO',
    'IIL':'ILE', 'IYR':'TYR', 'KCX':'LYS', 'LLP':'LYS', 'LLY':'LYS', 'LTR':'TRP', 'LYM':'LYS', 'LYZ':'LYS', 'MAA':'ALA', 'MEN':'ASN',
    'MHS':'HIS', 'MIS':'SER', 'MLE':'LEU', 'MPQ':'GLY', 'MSA':'GLY', 'MSE':'MET', 'MVA':'VAL', 'NEM':'HIS', 'NEP':'HIS', 'NLE':'LEU',
    'NLN':'LEU', 'NLP':'LEU', 'NMC':'GLY', 'OAS':'SER', 'OCS':'CYS', 'OMT':'MET', 'PAQ':'TYR', 'PCA':'GLU', 'PEC':'CYS', 'PHI':'PHE',
    'PHL':'PHE', 'PR3':'CYS', 'PRR':'ALA', 'PTR':'TYR', 'PYX':'CYS', 'SAC':'SER', 'SAR':'GLY', 'SCH':'CYS', 'SCS':'CYS', 'SCY':'CYS',
    'SEL':'SER', 'SEP':'SER', 'SET':'SER', 'SHC':'CYS', 'SHR':'LYS', 'SMC':'CYS', 'SOC':'CYS', 'STY':'TYR', 'SVA':'SER', 'TIH':'ALA',
    'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TRG':'LYS', 'TRO':'TRP', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR'
}

def _overlayPoints(points1, points2):
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

class PDBFixer(object):
    def __init__(self, pdb):
        self.pdb = pdb
        self.topology = pdb.topology
        self.positions = pdb.positions
        
        # Load the templates.
        
        self.templates = {}
        templatesPath = os.path.join(os.path.dirname(__file__), 'templates')
        for file in os.listdir(templatesPath):
            templatePdb = app.PDBFile(os.path.join(templatesPath, file))
            name = templatePdb.topology.residues().next().name
            self.templates[name] = templatePdb
        
    def _addAtomsToTopology(self, missingAtoms, missingTerminals, heavyAtomsOnly, omitUnknownMolecules):
        """Create a new Topology in which missing atoms have been added."""
        
        newTopology = app.Topology()
        newPositions = []*unit.nanometer
        newAtoms = []
        existingAtomMap = {}
        addedAtomMap = {}
        addedOXT = []
        for chain in self.topology.chains():
            if omitUnknownMolecules and not any(residue.name in self.templates for residue in chain.residues()):
                continue
            newChain = newTopology.addChain()
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain)
                
                # Add the existing heavy atoms.
                
                for atom in residue.atoms():
                    if not heavyAtomsOnly or (atom.element is not None and atom.element != hydrogen):
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        existingAtomMap[atom] = newAtom
                        newPositions.append(self.positions[atom.index])
                if residue in missingAtoms:
                    
                    # Find corresponding atoms in the residue and the template.
                    
                    template = self.templates[residue.name]
                    atomPositions = dict((atom.name, self.positions[atom.index]) for atom in residue.atoms())
                    points1 = []
                    points2 = []
                    for atom in template.topology.atoms():
                        if atom.name in atomPositions:
                            points1.append(atomPositions[atom.name].value_in_unit(unit.nanometer))
                            points2.append(template.positions[atom.index].value_in_unit(unit.nanometer))
                    
                    # Compute the optimal transform to overlay them.
                    
                    (translate2, rotate, translate1) = _overlayPoints(points1, points2)
                    
                    # Add the missing atoms.
                    
                    addedAtomMap[residue] = {}
                    for atom in missingAtoms[residue]:
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        newAtoms.append(newAtom)
                        addedAtomMap[residue][atom] = newAtom
                        templatePosition = template.positions[atom.index].value_in_unit(unit.nanometer)
                        newPositions.append((mm.Vec3(*np.dot(rotate, templatePosition+translate2))+translate1)*unit.nanometer)
                
                # If a terminal OXT is missing, add it.
                
                if residue in missingTerminals:
                    atomPositions = dict((atom.name, newPositions[atom.index].value_in_unit(unit.nanometer)) for atom in newResidue.atoms())
                    if 'OXT' in missingTerminals[residue]:
                        newAtom = newTopology.addAtom('OXT', oxygen, newResidue)
                        newAtoms.append(newAtom)
                        addedOXT.append(newAtom)
                        d_ca_o = atomPositions['O']-atomPositions['CA']
                        d_ca_c = atomPositions['C']-atomPositions['CA']
                        d_ca_c /= unit.sqrt(unit.dot(d_ca_c, d_ca_c))
                        v = d_ca_o - d_ca_c*unit.dot(d_ca_c, d_ca_o)
                        newPositions.append((atomPositions['O']+2*v)*unit.nanometer)
        
        # Add bonds from the original Topology.
        
        for atom1, atom2 in self.topology.bonds():
            if atom1 in existingAtomMap and atom2 in existingAtomMap:
                newTopology.addBond(existingAtomMap[atom1], existingAtomMap[atom2])
        
        # Add bonds that connect to new atoms.
        
        for residue in missingAtoms:
            template = self.templates[residue.name]
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
        
        # Return the results.
        
        return (newTopology, newPositions, newAtoms, existingAtomMap)
    
    def findNonstandardResidues(self):
        return [r for r in self.topology.residues() if r.name in substitutions]
    
    def replaceNonstandardResidues(self, replaceResidues):
        if len(replaceResidues) > 0:
            deleteAtoms = []

            # Find atoms that should be deleted.
            
            for residue in replaceResidues:
                residue.name = substitutions[residue.name]
                template = self.templates[residue.name]
                standardAtoms = set(atom.name for atom in template.topology.atoms())
                for atom in residue.atoms():
                    if atom.element in (None, hydrogen) or atom.name not in standardAtoms:
                        deleteAtoms.append(atom)
            
            # Delete them.
            
            modeller = app.Modeller(self.topology, self.positions)
            modeller.delete(deleteAtoms)
            self.topology = modeller.topology
            self.positions = modeller.positions
    
    def findMissingAtoms(self):
        missingAtoms = {}
        missingTerminals = {}
        
        # Loop over residues.
        
        for chain in self.topology.chains():
            chainResidues = list(chain.residues())
            for residue in chain.residues():
                if residue.name in self.templates:
                    template = self.templates[residue.name]
                    atomNames = set(atom.name for atom in residue.atoms())
                    
                    # Add atoms from the template that are missing.
                    
                    missing = []
                    for atom in template.topology.atoms():
                        if atom.name not in atomNames:
                            missing.append(atom)
                    if len(missing) > 0:
                        missingAtoms[residue] = missing
                    
                    # Add missing terminal atoms.
                    
                    terminals = []
                    if residue == chainResidues[-1]:
                        templateNames = set(atom.name for atom in template.topology.atoms())
                        if 'OXT' not in atomNames and all(name in templateNames for name in ['C', 'O', 'CA']):
                            terminals.append('OXT')
                        if len(terminals) > 0:
                            missingTerminals[residue] = terminals
        return (missingAtoms, missingTerminals)
    
    def addMissingAtoms(self, missingAtoms, missingTerminals):
        # Create a Topology that 1) adds missing atoms, 2) removes all hydrogens, and 3) removes unknown molecules.
        
        (newTopology, newPositions, newAtoms, existingAtomMap) = self._addAtomsToTopology(missingAtoms, missingTerminals, True, True)
        if len(newAtoms) > 0:
            
            # Create a System for energy minimizing it.
            
            res = list(newTopology.residues())
            forcefield = self._createForceField(newTopology)
            system = forcefield.createSystem(newTopology)
            
            # Set any previously existing atoms to be massless, they so won't move.
            
            for atom in existingAtomMap.itervalues():
                system.setParticleMass(atom.index, 0.0)
            
            # If any heavy atoms were omitted, add them back to avoid steric clashes.
            
            nonbonded = [f for f in system.getForces() if isinstance(f, mm.CustomNonbondedForce)][0]
            for atom in self.topology.atoms():
                if atom.element not in (None, hydrogen) and atom not in existingAtomMap:
                    system.addParticle(0.0)
                    nonbonded.addParticle([])
                    newPositions.append(self.positions[atom.index])
            
            # For efficiency, only compute interactions that involve a new atom.
            
            nonbonded.addInteractionGroup([atom.index for atom in newAtoms], range(system.getNumParticles()))
            
            # Do an energy minimization.
            
            integrator = mm.LangevinIntegrator(300*unit.kelvin, 10/unit.picosecond, 5*unit.femtosecond)
            context = mm.Context(system, integrator)
            context.setPositions(newPositions)
            mm.LocalEnergyMinimizer.minimize(context)
            state = context.getState(getPositions=True)
            
            # Now create a new Topology, including all atoms from the original one and adding the missing atoms.
            
            (newTopology2, newPositions2, newAtoms2, existingAtomMap2) = self._addAtomsToTopology(missingAtoms, missingTerminals, False, False)
            
            # Copy over the minimized positions for the new atoms.
            
            for a1, a2 in zip(newAtoms, newAtoms2):
                newPositions2[a2.index] = state.getPositions()[a1.index]
            self.topology = newTopology2
            self.positions = newPositions2
    
    def _createForceField(self, newTopology):
        forcefield = app.ForceField(os.path.join(os.path.dirname(__file__), 'soft.xml'))
        
        # The Topology may contain residues for which the ForceField does not have a template.
        # If so, we need to create new templates for them.
        
        atomTypes = {}
        bondedToAtom = []
        for atom in newTopology.atoms():
            bondedToAtom.append(set())
        for atom1, atom2 in newTopology.bonds():
            bondedToAtom[atom1.index].add(atom2.index)
            bondedToAtom[atom2.index].add(atom1.index)
        for residue in newTopology.residues():
            
            # Make sure the ForceField has a template for this residue.
            
            signature = app.forcefield._createResidueSignature([atom.element for atom in residue.atoms()])
            if signature in forcefield._templateSignatures:
                if any(app.forcefield._matchResidue(residue, t, bondedToAtom) is not None for t in forcefield._templateSignatures[signature]):
                    continue
            
            # Create a new template.
            
            resName = "extra_"+residue.name
            template = app.ForceField._TemplateData(resName)
            forcefield._templates[resName] = template
            indexInResidue = {}
            for atom in residue.atoms():
                element = atom.element
                typeName = 'extra_'+element.symbol
                if element not in atomTypes:
                    atomTypes[element] = (typeName, 0.0, element)
                    forcefield._atomTypes[typeName] = atomTypes[element]
                indexInResidue[atom.index] = len(template.atoms)
                template.atoms.append(app.ForceField._TemplateAtomData(atom.name, typeName, element))
            for atom in residue.atoms():
                for bondedTo in bondedToAtom[atom.index]:
                    if bondedTo in indexInResidue:
                        b = (indexInResidue[atom.index], indexInResidue[bondedTo])
                        if b[0] < b[1]:
                            template.bonds.append(b)
                            template.atoms[b[0]].bondedTo.append(b[1])
                            template.atoms[b[1]].bondedTo.append(b[0])
                    else:
                        b = indexInResidue[atom.index]
                        template.externalBonds.append(b)
                        template.atoms[b].externalBonds += 1
            if signature in forcefield._templateSignatures:
                forcefield._templateSignatures[signature].append(template)
            else:
                forcefield._templateSignatures[signature] = [template]
        return forcefield


if __name__=='__main__':
    if len(sys.argv) < 2:
        ui.launchUI()
    else:
        fixer = PDBFixer(app.PDBFile(sys.argv[1]))
        fixer.replaceNonstandardResidues(fixer.findNonstandardResidues())
        missingAtoms, missingTerminals = fixer.findMissingAtoms()
        fixer.addMissingAtoms(missingAtoms, missingTerminals)
        app.PDBFile.writeFile(fixer.topology, fixer.positions, open('output.pdb', 'w'))
