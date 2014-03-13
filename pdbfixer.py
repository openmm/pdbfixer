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
from simtk.openmm.app.internal.pdbstructure import PdbStructure
from simtk.openmm.app.element import hydrogen, oxygen
from simtk.openmm.app.forcefield import NonbondedGenerator
import numpy as np
import numpy.linalg as lin
import sys
import os
import os.path
import math

substitutions = {
    '2AS':'ASP', '3AH':'HIS', '5HP':'GLU', 'ACL':'ARG', 'AGM':'ARG', 'AIB':'ALA', 'ALM':'ALA', 'ALO':'THR', 'ALY':'LYS', 'ARM':'ARG',
    'ASA':'ASP', 'ASB':'ASP', 'ASK':'ASP', 'ASL':'ASP', 'ASQ':'ASP', 'AYA':'ALA', 'BCS':'CYS', 'BHD':'ASP', 'BMT':'THR', 'BNN':'ALA',
    'BUC':'CYS', 'BUG':'LEU', 'C5C':'CYS', 'C6C':'CYS', 'CAS':'CYS', 'CCS':'CYS', 'CEA':'CYS', 'CGU':'GLU', 'CHG':'ALA', 'CLE':'LEU', 'CME':'CYS',
    'CSD':'ALA', 'CSO':'CYS', 'CSP':'CYS', 'CSS':'CYS', 'CSW':'CYS', 'CSX':'CYS', 'CXM':'MET', 'CY1':'CYS', 'CY3':'CYS', 'CYG':'CYS',
    'CYM':'CYS', 'CYQ':'CYS', 'DAH':'PHE', 'DAL':'ALA', 'DAR':'ARG', 'DAS':'ASP', 'DCY':'CYS', 'DGL':'GLU', 'DGN':'GLN', 'DHA':'ALA',
    'DHI':'HIS', 'DIL':'ILE', 'DIV':'VAL', 'DLE':'LEU', 'DLY':'LYS', 'DNP':'ALA', 'DPN':'PHE', 'DPR':'PRO', 'DSN':'SER', 'DSP':'ASP',
    'DTH':'THR', 'DTR':'TRP', 'DTY':'TYR', 'DVA':'VAL', 'EFC':'CYS', 'FLA':'ALA', 'FME':'MET', 'GGL':'GLU', 'GL3':'GLY', 'GLZ':'GLY',
    'GMA':'GLU', 'GSC':'GLY', 'HAC':'ALA', 'HAR':'ARG', 'HIC':'HIS', 'HIP':'HIS', 'HMR':'ARG', 'HPQ':'PHE', 'HTR':'TRP', 'HYP':'PRO',
    'IAS':'ASP', 'IIL':'ILE', 'IYR':'TYR', 'KCX':'LYS', 'LLP':'LYS', 'LLY':'LYS', 'LTR':'TRP', 'LYM':'LYS', 'LYZ':'LYS', 'MAA':'ALA', 'MEN':'ASN',
    'MHS':'HIS', 'MIS':'SER', 'MLE':'LEU', 'MPQ':'GLY', 'MSA':'GLY', 'MSE':'MET', 'MVA':'VAL', 'NEM':'HIS', 'NEP':'HIS', 'NLE':'LEU',
    'NLN':'LEU', 'NLP':'LEU', 'NMC':'GLY', 'OAS':'SER', 'OCS':'CYS', 'OMT':'MET', 'PAQ':'TYR', 'PCA':'GLU', 'PEC':'CYS', 'PHI':'PHE',
    'PHL':'PHE', 'PR3':'CYS', 'PRR':'ALA', 'PTR':'TYR', 'PYX':'CYS', 'SAC':'SER', 'SAR':'GLY', 'SCH':'CYS', 'SCS':'CYS', 'SCY':'CYS',
    'SEL':'SER', 'SEP':'SER', 'SET':'SER', 'SHC':'CYS', 'SHR':'LYS', 'SMC':'CYS', 'SOC':'CYS', 'STY':'TYR', 'SVA':'SER', 'TIH':'ALA',
    'TPL':'TRP', 'TPO':'THR', 'TPQ':'ALA', 'TRG':'LYS', 'TRO':'TRP', 'TYB':'TYR', 'TYI':'TYR', 'TYQ':'TYR', 'TYS':'TYR', 'TYY':'TYR'
}
proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
rnaResidues = ['A', 'G', 'C', 'U', 'I']
dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']

def _overlayPoints(points1, points2):
    """Given two sets of points, determine the translation and rotation that matches them as closely as possible.
    
    This is based on W. Kabsch, Acta Cryst., A34, pp. 828-829 (1978)."""
    
    if len(points1) == 0:
        return (Vec3(0, 0, 0), np.identity(3), Vec3(0, 0, 0))
    if len(points1) == 1:
        return (points1[0], np.identity(3), -1*points2[0])
    
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
    """PDBFixer implements many tools for fixing problems in PDB files."""
    
    def __init__(self, structure):
        """Create a new PDBFixer to fix problems in a PDB file.
        
        Parameters:
         - structure (PdbStructure) the starting PDB structure containing problems to be fixed
        """
        self.structure = structure
        self.pdb = app.PDBFile(structure)
        self.topology = self.pdb.topology
        self.positions = self.pdb.positions
        self.centroid = unit.sum(self.positions)/len(self.positions)
        self.structureChains = list(self.structure.iter_chains())
        
        # Load the templates.
        
        self.templates = {}
        templatesPath = os.path.join(os.path.dirname(__file__), 'templates')
        for file in os.listdir(templatesPath):
            templatePdb = app.PDBFile(os.path.join(templatesPath, file))
            name = next(templatePdb.topology.residues()).name
            self.templates[name] = templatePdb
        
    def _addAtomsToTopology(self, heavyAtomsOnly, omitUnknownMolecules):
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
            chainResidues = list(chain.residues())
            newChain = newTopology.addChain()
            for indexInChain, residue in enumerate(chain.residues()):
                
                # Insert missing residues here.
                
                if (chain.index, indexInChain) in self.missingResidues:
                    insertHere = self.missingResidues[(chain.index, indexInChain)]
                    endPosition = self._computeResidueCenter(residue)
                    if indexInChain > 0:
                        startPosition = self._computeResidueCenter(chainResidues[indexInChain-1])
                    else:
                        outward = endPosition-self.centroid
                        norm = unit.norm(outward)
                        if norm > 0*unit.nanometer:
                            outward *= len(insertHere)*0.5*unit.nanometer/norm
                        startPosition = endPosition+outward
                    self._addMissingResiduesToChain(newChain, insertHere, startPosition, endPosition, residue, newAtoms, newPositions)
                
                # Create the new residue and add existing heavy atoms.
                                
                newResidue = newTopology.addResidue(residue.name, newChain)
                addResiduesAfter = (residue == chainResidues[-1] and (chain.index, indexInChain+1) in self.missingResidues)
                for atom in residue.atoms():
                    if not heavyAtomsOnly or (atom.element is not None and atom.element != hydrogen):
                        if atom.name == 'OXT' and (chain.index, indexInChain+1) in self.missingResidues:
                            continue # Remove terminal oxygen, since we'll add more residues after this one
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        existingAtomMap[atom] = newAtom
                        newPositions.append(self.positions[atom.index])
                if residue in self.missingAtoms:
                    
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
                    for atom in self.missingAtoms[residue]:
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        newAtoms.append(newAtom)
                        addedAtomMap[residue][atom] = newAtom
                        templatePosition = template.positions[atom.index].value_in_unit(unit.nanometer)
                        newPositions.append((mm.Vec3(*np.dot(rotate, templatePosition+translate2))+translate1)*unit.nanometer)
                if residue in self.missingTerminals:
                    terminalsToAdd = self.missingTerminals[residue]
                else:
                    terminalsToAdd = None

                # If this is the end of the chain, add any missing residues that come after it.
                
                if residue == chainResidues[-1] and (chain.index, indexInChain+1) in self.missingResidues:
                    insertHere = self.missingResidues[(chain.index, indexInChain+1)]
                    if len(insertHere) > 0:
                        startPosition = self._computeResidueCenter(residue)
                        outward = startPosition-self.centroid
                        norm = unit.norm(outward)
                        if norm > 0*unit.nanometer:
                            outward *= len(insertHere)*0.5*unit.nanometer/norm
                        endPosition = startPosition+outward
                        self._addMissingResiduesToChain(newChain, insertHere, startPosition, endPosition, residue, newAtoms, newPositions)
                        newResidue = list(newChain.residues())[-1]
                        if newResidue.name in proteinResidues:
                            terminalsToAdd = ['OXT']
                        else:
                            terminalsToAdd = None
                
                # If a terminal OXT is missing, add it.
                
                if terminalsToAdd is not None:
                    atomPositions = dict((atom.name, newPositions[atom.index].value_in_unit(unit.nanometer)) for atom in newResidue.atoms())
                    if 'OXT' in terminalsToAdd:
                        newAtom = newTopology.addAtom('OXT', oxygen, newResidue)
                        newAtoms.append(newAtom)
                        addedOXT.append(newAtom)
                        d_ca_o = atomPositions['O']-atomPositions['CA']
                        d_ca_c = atomPositions['C']-atomPositions['CA']
                        d_ca_c /= unit.sqrt(unit.dot(d_ca_c, d_ca_c))
                        v = d_ca_o - d_ca_c*unit.dot(d_ca_c, d_ca_o)
                        newPositions.append((atomPositions['O']+2*v)*unit.nanometer)
        newTopology.setUnitCellDimensions(self.topology.getUnitCellDimensions())
        newTopology.createStandardBonds()
        newTopology.createDisulfideBonds(newPositions)
        
        # Return the results.
        
        return (newTopology, newPositions, newAtoms, existingAtomMap)
    
    def _computeResidueCenter(self, residue):
        """Compute the centroid of a residue."""
        return unit.sum([self.pdb.positions[atom.index] for atom in residue.atoms()])/len(list(residue.atoms()))
    
    def _addMissingResiduesToChain(self, chain, residueNames, startPosition, endPosition, orientTo, newAtoms, newPositions):
        """Add a series of residues to a chain."""
        orientToPositions = dict((atom.name, self.positions[atom.index]) for atom in orientTo.atoms())
        for i, residueName in enumerate(residueNames):
            template = self.templates[residueName]

            # Find a translation that best matches the adjacent residue.
        
            points1 = []
            points2 = []
            for atom in template.topology.atoms():
                if atom.name in orientToPositions:
                    points1.append(orientToPositions[atom.name].value_in_unit(unit.nanometer))
                    points2.append(template.positions[atom.index].value_in_unit(unit.nanometer))
            (translate2, rotate, translate1) = _overlayPoints(points1, points2)
            
            # Create the new residue.
            
            newResidue = chain.topology.addResidue(residueName, chain)
            translate = startPosition+(endPosition-startPosition)*(i+1.0)/(len(residueNames)+1.0)
            templateAtoms = list(template.topology.atoms())
            if newResidue == next(chain.residues()):
                templateAtoms = [atom for atom in templateAtoms if atom.name not in ('P', 'OP1', 'OP2')]
            for atom in templateAtoms:
                newAtom = chain.topology.addAtom(atom.name, atom.element, newResidue)
                newAtoms.append(newAtom)
                templatePosition = template.positions[atom.index].value_in_unit(unit.nanometer)
                newPositions.append(mm.Vec3(*np.dot(rotate, templatePosition))*unit.nanometer+translate)
    
    def removeChains(self, chainIndices):
        """Remove a set of chains from the structure.
        
        Parameters:
         - chainIndices (list) the indices of the chains to remove.
        """
        modeller = app.Modeller(self.topology, self.positions)
        allChains = list(self.topology.chains())
        modeller.delete(allChains[i] for i in chainIndices)
        self.topology = modeller.topology
        self.positions = modeller.positions
        self.structureChains = [self.structureChains[i] for i in range(len(self.structureChains)) if i not in chainIndices]
    
    def findMissingResidues(self):
        """Find residues that are missing from the structure.
        
        The results are stored into the missingResidues field, which is a dict.  Each key is a tuple consisting of
        the index of a chain, and the residue index within that chain at which new residues should be inserted.
        The corresponding value is a list of the names of residues to insert there.
        """
        chains = [c for c in self.structureChains if any(atom.record_name == 'ATOM' for atom in c.iter_atoms())]
        chainWithGaps = {}
        
        # Find the sequence of each chain, with gaps for missing residues.
        
        for chain in chains:
            minResidue = min(r.number for r in chain.iter_residues())
            maxResidue = max(r.number for r in chain.iter_residues())
            residues = [None]*(maxResidue-minResidue+1)
            for r in chain.iter_residues():
                residues[r.number-minResidue] = r.get_name()
            chainWithGaps[chain] = residues
        
        # Try to find the chain that matches each sequence.
        
        chainSequence = {}
        chainOffset = {}
        for sequence in self.structure.sequences:
            for chain in chains:
                if chain.chain_id != sequence.chain_id:
                    continue
                if chain in chainSequence:
                    continue
                for offset in range(len(sequence.residues)-len(chainWithGaps[chain])+1):
                    if all(a == b or b == None for a,b in zip(sequence.residues[offset:], chainWithGaps[chain])):
                        chainSequence[chain] = sequence
                        chainOffset[chain] = offset
                        break
                if chain in chainSequence:
                    break
        
        # Now build the list of residues to add.
        
        self.missingResidues = {}
        for structChain, topChain in zip(self.structureChains, self.topology.chains()):
            if structChain in chainSequence:
                offset = chainOffset[structChain]
                sequence = chainSequence[structChain].residues
                gappedSequence = chainWithGaps[structChain]
                index = 0
                for i in range(len(sequence)):
                    if i < offset or i >= len(gappedSequence)+offset or gappedSequence[i-offset] is None:
                        key = (topChain.index, index)
                        if key not in self.missingResidues:
                            self.missingResidues[key] = []
                        residueName = sequence[i]
                        if residueName in substitutions:
                            residueName = substitutions[sequence[i]]
                        self.missingResidues[key].append(residueName)
                    else:
                        index += 1
    
    def findNonstandardResidues(self):
        """Identify non-standard residues found in the structure, and select standard residues to replace them with.
        
        The results are stored into the nonstandardResidues field, which is a map of Residue objects to the names
        of suggested replacement residues.
        """
        
        # First find residues based on our table of standard substitutions.
        
        nonstandard = dict((r, substitutions[r.name]) for r in self.topology.residues() if r.name in substitutions)
        
        # Now add ones based on MODRES records.
        
        modres = dict(((m.chain_id, m.number, m.residue_name), m.standard_name) for m in self.structure.modified_residues)
        for structChain, topChain in zip(self.structureChains, self.topology.chains()):
            for structResidue, topResidue in zip(structChain.iter_residues(), topChain.residues()):
                key = (structChain.chain_id, structResidue.number, structResidue.name)
                if key in modres:
                    replacement = modres[key]
                    if replacement == 'DU':
                        replacement = 'DT'
                    if replacement in self.templates:
                        nonstandard[topResidue] = replacement
        self.nonstandardResidues = [(r, nonstandard[r]) for r in sorted(nonstandard, key=lambda r: r.index)]
    
    def replaceNonstandardResidues(self):
        """Replace every residue listed in the nonstandardResidues field with the specified standard residue."""
        if len(self.nonstandardResidues) > 0:
            deleteAtoms = []

            # Find atoms that should be deleted.
            
            for residue, replaceWith in self.nonstandardResidues:
                residue.name = replaceWith
                template = self.templates[replaceWith]
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
        """Find heavy atoms that are missing from the structure.
        
        The results are stored into two fields: missingAtoms and missingTerminals.  Each of these is a dict whose keys
        are Residue objects and whose values are lists of atom names.  missingAtoms contains standard atoms that should
        be present in any residue of that type.  missingTerminals contains terminal atoms that should be present at the
        start or end of a chain.
        """
        missingAtoms = {}
        missingTerminals = {}
        
        # Loop over residues.
        
        for chain in self.topology.chains():
            chainResidues = list(chain.residues())
            for residue in chain.residues():
                if residue.name in self.templates:
                    template = self.templates[residue.name]
                    atomNames = set(atom.name for atom in residue.atoms())
                    templateAtoms = list(template.topology.atoms())
                    if residue == chainResidues[0] and (chain.index, 0) not in self.missingResidues:
                        templateAtoms = [atom for atom in templateAtoms if atom.name not in ('P', 'OP1', 'OP2')]
                    
                    # Add atoms from the template that are missing.
                    
                    missing = []
                    for atom in templateAtoms:
                        if atom.name not in atomNames:
                            missing.append(atom)
                    if len(missing) > 0:
                        missingAtoms[residue] = missing
                    
                    # Add missing terminal atoms.
                    
                    terminals = []
                    if residue == chainResidues[-1] and (chain.index, len(chainResidues)) not in self.missingResidues:
                        templateNames = set(atom.name for atom in template.topology.atoms())
                        if 'OXT' not in atomNames and all(name in templateNames for name in ['C', 'O', 'CA']):
                            terminals.append('OXT')
                        if len(terminals) > 0:
                            missingTerminals[residue] = terminals
        self.missingAtoms = missingAtoms
        self.missingTerminals = missingTerminals
    
    def addMissingAtoms(self):
        """Add all missing heavy atoms, as specified by the missingAtoms, missingTerminals, and missingResidues fields."""
        
        # Create a Topology that 1) adds missing atoms, 2) removes all hydrogens, and 3) removes unknown molecules.
        
        (newTopology, newPositions, newAtoms, existingAtomMap) = self._addAtomsToTopology(True, True)
        if len(newAtoms) > 0:
            
            # Create a System for energy minimizing it.
            
            res = list(newTopology.residues())
            forcefield = self._createForceField(newTopology, False)
            system = forcefield.createSystem(newTopology)
            
            # Set any previously existing atoms to be massless, they so won't move.
            
            for atom in existingAtomMap.values():
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
            nearest = self._findNearestDistance(context, newTopology, newAtoms)
            if nearest < 0.15:
                
                # Some atoms are very close together.  Run some dynamics while slowly increasing the strength of the
                # repulsive interaction to try to improve the result.
                
                for i in range(10):
                    context.setParameter('C', 0.15*(i+1))
                    integrator.step(200)
                    d = self._findNearestDistance(context, newTopology, newAtoms)
                    if d > nearest:
                        nearest = d
                        state = context.getState(getPositions=True)
                        if nearest >= 0.15:
                            break
                context.setState(state)
                state = context.getState(getPositions=True)
            
            # Now create a new Topology, including all atoms from the original one and adding the missing atoms.
            
            (newTopology2, newPositions2, newAtoms2, existingAtomMap2) = self._addAtomsToTopology(False, False)
            
            # Copy over the minimized positions for the new atoms.
            
            for a1, a2 in zip(newAtoms, newAtoms2):
                newPositions2[a2.index] = state.getPositions()[a1.index]
            self.topology = newTopology2
            self.positions = newPositions2
    
    def removeHeterogens(self, keepWater):
        """Remove all heterogens from the structure.
        
        Parameters:
         - keepWater (bool) if True, water molecules will not be removed
        """
        keep = set(proteinResidues).union(dnaResidues).union(rnaResidues)
        keep.add('N')
        keep.add('UNK')
        if keepWater:
            keep.add('HOH')
        toDelete = []
        for residue in self.topology.residues():
            if residue.name not in keep:
                toDelete.append(residue)
        modeller = app.Modeller(self.topology, self.positions)
        modeller.delete(toDelete)
        self.topology = modeller.topology
        self.positions = modeller.positions
    
    def addMissingHydrogens(self, pH):
        """Add missing hydrogen atoms to the structure.
        
        Parameters:
         - pH (float) the pH based on which to select hydrogens
        """
        modeller = app.Modeller(self.topology, self.positions)
        modeller.addHydrogens(pH=pH)
        self.topology = modeller.topology
        self.positions = modeller.positions
    
    def addSolvent(self, boxSize, positiveIon='Na+', negativeIon='Cl-', ionicStrength=0*unit.molar):
        """Add a solvent box surrounding the structure.
        
        Parameters:
         - boxSize (Vec3) the size of the box to fill with water
         - positiveIon (string='Na+') the type of positive ion to add.  Allowed values are 'Cs+', 'K+', 'Li+', 'Na+', and 'Rb+'
         - negativeIon (string='Cl-') the type of negative ion to add.  Allowed values are 'Cl-', 'Br-', 'F-', and 'I-'
         - ionicString (concentration=0*molar) the total concentration of ions (both positive and negative) to add.  This
           does not include ions that are added to neutralize the system.
        """
        modeller = app.Modeller(self.topology, self.positions)
        forcefield = self._createForceField(self.topology, True)
        modeller.addSolvent(forcefield, boxSize=boxSize, positiveIon=positiveIon, negativeIon=negativeIon, ionicStrength=ionicStrength)
        self.topology = modeller.topology
        self.positions = modeller.positions
    
    def _createForceField(self, newTopology, water):
        """Create a force field to use for optimizing the positions of newly added atoms."""
        
        if water:
            forcefield = app.ForceField('amber10.xml', 'tip3p.xml')
            nonbonded = [f for f in forcefield._forces if isinstance(f, NonbondedGenerator)][0]
            radii = {'H':0.198, 'Li':0.203, 'C':0.340, 'N':0.325, 'O':0.299, 'F':0.312, 'Na':0.333, 'Mg':0.141,
                     'P':0.374, 'S':0.356, 'Cl':0.347, 'K':0.474, 'Br':0.396, 'Rb':0.527, 'I':0.419, 'Cs':0.605}
        else:
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
                    if water:
                        # Select a reasonable vdW radius for this atom type.
                        
                        if element.symbol in radii:
                            sigma = radii[element.symbol]
                        else:
                            sigma = 0.5
                        nonbonded.typeMap[typeName] = (0.0, sigma, 0.0)
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
    
    def _findNearestDistance(self, context, topology, newAtoms):
        """Given a set of newly added atoms, find the closest distance between one of those atoms and another atom."""
        
        positions = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        atomResidue = [atom.residue for atom in topology.atoms()]
        nearest = sys.float_info.max
        for atom in newAtoms:
            p = positions-positions[atom.index]
            dist = math.sqrt(min(np.dot(p[i], p[i]) for i in range(len(atomResidue)) if atomResidue[i] != atom.residue))
            if dist < nearest:
                nearest = dist
        return nearest


if __name__=='__main__':
    if len(sys.argv) < 2:
        # Display the UI.
        
        import ui
        ui.launchUI()
    else:
        # Run in command line mode.
        
        from optparse import OptionParser
        parser = OptionParser(usage="Usage: %prog\n       %prog [options] filename\n\nWhen run with no arguments, it launches the user interface.  If any arguments are specified, it runs in command line mode.")
        parser.add_option('--output', default='output.pdb', dest='output', metavar='FILENAME', help='output pdb file [default: output.pdb]')
        parser.add_option('--add-atoms', default='all', dest='atoms', choices=('all', 'heavy', 'hydrogen', 'none'), help='which missing atoms to add: all, heavy, hydrogen, or none [default: all]')
        parser.add_option('--keep-heterogens', default='all', dest='heterogens', choices=('all', 'water', 'none'), metavar='OPTION', help='which heterogens to keep: all, water, or none [default: all]')
        parser.add_option('--replace-nonstandard', action='store_true', default=False, dest='nonstandard', help='replace nonstandard residues with standard equivalents')
        parser.add_option('--add-residues', action='store_true', default=False, dest='residues', help='add missing residues')
        parser.add_option('--water-box', dest='box', type='float', nargs=3, metavar='X Y Z', help='add a water box. The value is the box dimensions in nm [example: --water-box=2.5 2.4 3.0]')
        parser.add_option('--ph', type='float', default=7.0, dest='ph', help='the pH to use for adding missing hydrogens [default: 7.0]')
        parser.add_option('--positive-ion', default='Na+', dest='positiveIon', choices=('Cs+', 'K+', 'Li+', 'Na+', 'Rb+'), metavar='ION', help='positive ion to include in the water box: Cs+, K+, Li+, Na+, or Rb+ [default: Na+]')
        parser.add_option('--negative-ion', default='Cl-', dest='negativeIon', choices=('Cl-', 'Br-', 'F-', 'I-'), metavar='ION', help='negative ion to include in the water box: Cl-, Br-, F-, or I- [default: Cl-]')
        parser.add_option('--ionic-strength', type='float', default=0.0, dest='ionic', metavar='STRENGTH', help='molar concentration of ions to add to the water box [default: 0.0]')
        (options, args) = parser.parse_args()
        if len(args) == 0:
            parser.error('No filename specified')
        if len(args) > 1:
            parser.error('Must specify a single filename')
        fixer = PDBFixer(PdbStructure(open(args[0])))
        if options.residues:
            fixer.findMissingResidues()
        else:
            fixer.missingResidues = {}
        if options.nonstandard:
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        if options.atoms not in ('all', 'heavy'):
            fixer.missingAtoms = {}
            fixer.missingTerminals = {}
        fixer.addMissingAtoms()
        if options.heterogens == 'none':
            fixer.removeHeterogens(False)
        elif options.heterogens == 'water':
            fixer.removeHeterogens(True)
        if options.atoms in ('all', 'hydrogen'):
            fixer.addMissingHydrogens(options.ph)
        if options.box is not None:
            fixer.addSolvent(options.box*unit.nanometer, options.positiveIon, options.negativeIon, options.ionic*unit.molar)
        app.PDBFile.writeFile(fixer.topology, fixer.positions, open(options.output, 'w'))
