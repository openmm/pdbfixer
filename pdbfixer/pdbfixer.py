"""
pdbfixer.py: Fixes problems in PDB files

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2013-2014 Stanford University and the Authors.
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
__version__ = "1.1"

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

# Imports for urlopen
if sys.version_info >= (3,0):
    from urllib.request import urlopen
else:
    from urllib2 import urlopen
        
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

    Parameters
    ----------
    points1 (numpy array of simtk.unit.Quantity with units compatible with distance) - reference set of coordinates
    points2 (numpy array of simtk.unit.Quantity with units compatible with distance) - set of coordinates to be rotated

    Returns
    -------
    translate2 - vector to translate points2 by in order to center it
    rotate - rotation matrix to apply to centered points2 to map it on to points1
    center1 - center of points1

    Notes
    -----    
    This is based on W. Kabsch, Acta Cryst., A34, pp. 828-829 (1978).

    """
    
    if len(points1) == 0:
        return (mm.Vec3(0, 0, 0), np.identity(3), mm.Vec3(0, 0, 0))
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

def _findUnoccupiedDirection(point, positions):
    """Given a point in space and a list of atom positions, find the direction in which the local density of atoms is lowest."""
    
    point = point.value_in_unit(unit.nanometers)
    direction = mm.Vec3(0, 0, 0)
    for pos in positions.value_in_unit(unit.nanometers):
        delta = pos-point
        distance = unit.norm(delta)
        if distance > 0.1:
            distance2 = distance*distance
            direction -= delta/(distance2*distance2)
    direction /= unit.norm(direction)
    return direction

class PDBFixer(object):
    """PDBFixer implements many tools for fixing problems in PDB files.
    """
    
    def __init__(self, filename=None, pdbfile=None, url=None, pdbid=None):
        """Create a new PDBFixer instance to fix problems in a PDB file.
        
        Parameters
        ----------
        filename : str, optional, default=None
            A filename specifying the file from which the PDB file is to be read.
        pdbfile : file, optional, default=None
            A file-like object from which the PDB file is to be read.
            The file is not closed after reading.
        url : str, optional, default=None
            A URL specifying the internet location from which the PDB file contents should be retrieved.
        pdbid : str, optional, default=None
            A four-letter PDB code specifying the structure to be retrieved from the RCSB.
            
        Notes
        -----
        Only one of structure, filename, pdbfile, url, or pdbid may be specified or an exception will be thrown.
            
        Examples
        --------
        
        Start from a file object.

        >>> pdbid = '1VII'
        >>> url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
        >>> file = urlopen(url)
        >>> fixer = PDBFixer(pdbfile=file)

        Start from a filename.
        
        >>> filename = 'test.pdb'
        >>> file = urlopen(url)
        >>> outfile = open(filename, 'w')
        >>> outfile.write(file.read())
        >>> outfile.close()
        >>> fixer = PDBFixer(filename=filename)
        
        Start from a URL.

        >>> fixer = PDBFixer(url=url)

        Start from a PDB code.
        
        >>> fixer = PDBFixer(pdbid=pdbid)

        """

        # Check to make sure only one option has been specified.
        if bool(filename) + bool(pdbfile) + bool(url) + bool(pdbid) != 1:
            raise Exception("Exactly one option [filename, pdbfile, url, pdbid] must be specified.")

        self.source = None
        if filename:
            self.source = filename
            # A local file has been specified.
            file = open(filename, 'r')                
            structure = PdbStructure(file)
            file.close()
        elif pdbfile:
            # A file-like object has been specified.
            structure = PdbStructure(pdbfile)  
        elif url:
            self.source = url
            # A URL has been specified.
            file = urlopen(url)
            structure = PdbStructure(file)
            file.close()
        elif pdbid:
            # A PDB id has been specified.
            url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbid
            self.source = url
            file = urlopen(url)
            # Read contents all at once and split into lines, since urlopen doesn't like it when we read one line at a time over the network.
            contents = file.read()
            lines = contents.split('\n')
            file.close()
            structure = PdbStructure(lines)
            
        # Check the structure has some atoms in it.
        atoms = list(structure.iter_atoms())
        if len(atoms)==0:
            raise Exception("Structure contains no atoms.")
            
        self.structure = structure
        self.pdb = app.PDBFile(structure)
        self.topology = self.pdb.topology
        self.positions = self.pdb.positions
        self.structureChains = list(self.structure.iter_chains())
        
        # Load the templates.
        
        self.templates = {}
        templatesPath = os.path.join(os.path.dirname(__file__), 'templates')
        for file in os.listdir(templatesPath):
            templatePdb = app.PDBFile(os.path.join(templatesPath, file))
            name = next(templatePdb.topology.residues()).name
            self.templates[name] = templatePdb
        
        return

    def _addAtomsToTopology(self, heavyAtomsOnly, omitUnknownMolecules):
        """Create a new Topology in which missing atoms have been added.

        Parameters
        ----------
        heavyAtomsOnly : bool
            If True, only heavy atoms will be added to the topology.
        omitUnknownMolecules : bool
            If True, unknown molecules will be omitted from the topology.

        Returns
        -------
        newTopology : simtk.openmm.app.Topology
            A new Topology object containing atoms from the old.         
        newPositions : list of simtk.unit.Quantity with units compatible with nanometers
            Atom positions for the new Topology object.
        newAtoms : simtk.openmm.app.Topology.Atom
            New atom objects.
        existingAtomMap : dict
            Mapping from old atoms to new atoms.

        """
        
        newTopology = app.Topology()
        newPositions = []*unit.nanometer
        newAtoms = []
        existingAtomMap = {}
        addedAtomMap = {}
        addedOXT = []
        residueCenters = [self._computeResidueCenter(res).value_in_unit(unit.nanometers) for res in self.topology.residues()]*unit.nanometers
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
                        loopDirection = _findUnoccupiedDirection((startPosition+endPosition)/2, residueCenters)
                    else:
                        outward = _findUnoccupiedDirection(endPosition, residueCenters)*unit.nanometers
                        norm = unit.norm(outward)
                        if norm > 0*unit.nanometer:
                            outward *= len(insertHere)*0.5*unit.nanometer/norm
                        startPosition = endPosition+outward
                        loopDirection = None
                    self._addMissingResiduesToChain(newChain, insertHere, startPosition, endPosition, loopDirection, residue, newAtoms, newPositions)
                
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
                        outward = _findUnoccupiedDirection(startPosition, residueCenters)*unit.nanometers
                        norm = unit.norm(outward)
                        if norm > 0*unit.nanometer:
                            outward *= len(insertHere)*0.5*unit.nanometer/norm
                        endPosition = startPosition+outward
                        self._addMissingResiduesToChain(newChain, insertHere, startPosition, endPosition, None, residue, newAtoms, newPositions)
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
    
    def _addMissingResiduesToChain(self, chain, residueNames, startPosition, endPosition, loopDirection, orientTo, newAtoms, newPositions):
        """Add a series of residues to a chain."""
        orientToPositions = dict((atom.name, self.positions[atom.index]) for atom in orientTo.atoms())
        if loopDirection is None:
            loopDirection = mm.Vec3(0, 0, 0)
        
        # We'll add the residues in an arc connecting the endpoints.  Figure out the height of that arc.
        
        length = unit.norm(endPosition-startPosition)
        numResidues = len(residueNames)
        if length > numResidues*0.3*unit.nanometers:
            loopHeight = 0*unit.nanometers
        else:
            loopHeight = (numResidues*0.3*unit.nanometers-length)/2
        
        # Add the residues.
        
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
            fraction = (i+1.0)/(numResidues+1.0)
            translate = startPosition + (endPosition-startPosition)*fraction + loopHeight*math.sin(fraction*math.pi)*loopDirection
            templateAtoms = list(template.topology.atoms())
            if newResidue == next(chain.residues()):
                templateAtoms = [atom for atom in templateAtoms if atom.name not in ('P', 'OP1', 'OP2')]
            for atom in templateAtoms:
                newAtom = chain.topology.addAtom(atom.name, atom.element, newResidue)
                newAtoms.append(newAtom)
                templatePosition = template.positions[atom.index].value_in_unit(unit.nanometer)
                newPositions.append(mm.Vec3(*np.dot(rotate, templatePosition))*unit.nanometer+translate)
    
    def removeChains(self, chainIndices=None, chainIds=None):
        """Remove a set of chains from the structure.
        
        Parameters
        ----------
        chainIndices : list of int, optional, default=None
            List of indices of chains to remove.
        chainIds : list of str, optional, default=None
            List of chain ids of chains to remove.

        Examples
        --------

        Load a PDB file with two chains and eliminate the second chain.

        >>> fixer = PDBFixer(pdbid='4J7F')
        >>> fixer.removeChains(chainIndices=[1])

        Load a PDB file with two chains and eliminate chains named 'B' and 'D'.

        >>> fixer = PDBFixer(pdbid='4J7F')
        >>> fixer.removeChains(chainIds=['B','D'])

        """
        modeller = app.Modeller(self.topology, self.positions)
        allChains = list(self.topology.chains())

        if chainIndices == None:
            chainIndices = list()
        if chainIds != None:
            # Add all chains that match the selection to the list.
            chain_id_list = [c.chain_id for c in self.structureChains]
            for (chainNumber, chainId) in enumerate(chain_id_list):
                if chainId in chainIds:
                    chainIndices.append(chainNumber)
            # Ensure only unique entries remain.
            chainIndices = list(set(chainIndices))

        # Do nothing if no chains will be deleted.
        if len(chainIndices) == 0:
            return

        modeller.delete(allChains[i] for i in chainIndices)
        self.topology = modeller.topology
        self.positions = modeller.positions
        self.structureChains = [self.structureChains[i] for i in range(len(self.structureChains)) if i not in chainIndices]

        return

    def removeResiduesFromChains(self, retainStart, retainEnd, chainIndices=None, chainIds=None):
        """Remove all residues except the specified span (inclusive).

        If no chains are specified, only the residue span is retained from all chains.

        Parameters
        ----------
        retainStart, retainEnd : int
            Initial and final residues (inclusive) in span to retain; all other residues from specified chains are removed.
        chainIndices : list of int, optional, default=None
            List of indices of chains from which residues are to be removed.
        chainIds : list of str, optional, default=None
            List of chain ids of chains from which residues are to be removed.

        Examples
        --------

        Load a PDB file with two chains and eliminate specified residues from both chains.

        >>> fixer = PDBFixer(pdbid='4JSV')
        >>> fixer.removeResiduesFromChains(2001, 2549)

        Load a PDB file with two chains and eliminate specified residues from chain B.

        >>> fixer = PDBFixer(pdbid='4JSV')
        >>> fixer.removeResiduesFromChains(2001, 2549, chainIds=['B'])

        Load a PDB file with two chains and eliminate specified residues from the second chain.

        >>> fixer = PDBFixer(pdbid='4JSV')
        >>> fixer.removeResiduesFromChains(2001, 2549, chainIndices=[1])

        """
        # Iterate of residues from structure and topology at the same time.
        structure_residues = list()
        for chain in self.structureChains:
            structure_residues += chain.residues
        topology_residues = [ residue for residue in self.topology.residues() ]

        # Sanity check: Make sure residue info matches up.
        for (structure_residue, topology_residue) in zip(structure_residues, topology_residues):
            if (structure_residue.name != topology_residue.name):
                raise Exception("structure_residue.name != topology_residue.name: structure_residue = [%d %s], topology_residue = [%d %s]" % (structure_residue.number, structure_residue.name, topology_residue.index, topology_residue.name))

        # Create mapping.
        structure_to_topology = { structure_residue : topology_residue for (structure_residue, topology_residue) in zip(structure_residues, topology_residues) }

        # Determine chains to process.
        if chainIndices == None:
            chainIndices = list()
        if chainIds != None:
            # Add all chains that match the selection to the list.
            chain_id_list = [c.chain_id for c in self.structureChains]
            for (chainNumber, chainId) in enumerate(chain_id_list):
                if chainId in chainIds:
                    chainIndices.append(chainNumber)
            # Ensure only unique entries remain.
            chainIndices = list(set(chainIndices))
        if (chainIndices==None) and (chainIds==None):
            # Add all chains.
            chain_id_list = [c.chain_id for c in self.structureChains]
            for (chainNumber, chainId) in enumerate(chain_id_list):
                chainIndices.append(chainNumber)

        # Build residue deletion queue.
        topology_residues_to_delete = list()
        structure_residues_to_delete = list()
        chains = [ c for c in self.structure.iter_chains() ]
        chains = [ chains[index] for index in chainIndices ]
        for chain in chains:
            for structure_residue in chain.residues:
                topology_residue = structure_to_topology[structure_residue]
                if (structure_residue.number < retainStart) or (structure_residue.number > retainEnd):
                    # Add to queue for deletion.
                    topology_residues_to_delete.append(topology_residue)
                    structure_residues_to_delete.append((chain, structure_residue))

        # Delete residues.
        modeller = app.Modeller(self.topology, self.positions)
        modeller.delete(topology_residues_to_delete)
        [self.topology, self.positions] = [modeller.topology, modeller.positions]

        # Delete structure residues.
        for (chain, structure_residue) in structure_residues_to_delete:
            chain.residues.remove(structure_residue)
            del chain.residues_by_num_icode[str(structure_residue.number) + structure_residue.insertion_code]
            del chain.residues_by_number[structure_residue.number]

        return

    def findMissingResidues(self):
        """Find residues that are missing from the structure.
        
        The results are stored into the missingResidues field, which is a dict.  Each key is a tuple consisting of
        the index of a chain, and the residue index within that chain at which new residues should be inserted.
        The corresponding value is a list of the names of residues to insert there.

        Examples
        --------

        >>> fixer = PDBFixer(pdbid='1VII')
        >>> fixer.findMissingResidues()
        >>> missing_residues = fixer.missingResidues

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

        Examples
        --------

        Find nonstandard residues.

        >>> fixer = PDBFixer(pdbid='1YRI')
        >>> fixer.findNonstandardResidues()
        >>> nonstandard_residues = fixer.nonstandardResidues        

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
        """Replace every residue listed in the nonstandardResidues field with the specified standard residue.

        Notes
        -----
        You must have first called findNonstandardResidues() to identify nonstandard residues.

        Examples
        --------

        Find and replace nonstandard residues using replacement templates stored in the 'templates' field of PDBFixer object.

        >>> fixer = PDBFixer(pdbid='1YRI')
        >>> fixer.findNonstandardResidues()
        >>> fixer.replaceNonstandardResidues()

        """
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


    def applyMutations(self, mutations, chain_id, which_model=0):
        """Apply a list of amino acid substitutions to make a mutant protein.

        Parameters
        ----------
        mutations : list of strings
            Each string must include the resName (original), index,
            and resName (target).  For example, ALA-133-GLY will mutate
            alanine 133 to glycine.
        chain_id : str
            String based chain ID of the single chain you wish to mutate.
        which_model : int, default = 0
            Which model to use in the pdb structure.

        Notes
        -----

        We require three letter codes to avoid possible ambiguitities.
        We can't guarnatee that the resulting model is a good one; for
        significant changes in sequence, you should probably be using
        a standalone homology modelling tool.

        Examples
        --------

        Apply a single mutation to chain A.

        >>> fixer = PDBFixer(pdbid='1VII')
        >>> fixer.applyMutations(["ALA-57-GLY"], chain_id='A')

        """
        # First find residues based on our table of standard substitutions.
        
        index_to_old_name = dict((r.index, r.name) for r in self.topology.residues())
        index_to_new_residues = {}
        
        chain_id_to_chain_number = dict((c.chain_id, k) for k, c in enumerate(self.structure.models[which_model].chains))
        chain_number = chain_id_to_chain_number[chain_id]
        
        resSeq_to_index = dict((r.number, k) for k, r in enumerate(self.structure.models[which_model].chains[chain_number]))
        
        for mut_str in mutations:
            old_name, resSeq, new_name = mut_str.split("-")
            resSeq = int(resSeq)
            index = resSeq_to_index[resSeq]
            
            if index not in index_to_old_name:
                raise(KeyError("Cannot find index %d in system!" % index))
            
            if index_to_old_name[index] != old_name:
                raise(ValueError("You asked to mutate %s %d, but that residue is actually %s!" % (old_name, index, index_to_old_name[index])))
            
            try:
                template = self.templates[new_name]
            except KeyError:
                raise(KeyError("Cannot find residue %s in template library!" % new_name))
            
            index_to_new_residues[index] = new_name
            
        
        residue_map = [(r, index_to_new_residues[r.index]) for r in self.topology.residues() if r.index in index_to_new_residues]

        if len(residue_map) > 0:
            deleteAtoms = []

            # Find atoms that should be deleted.
            
            for residue, replaceWith in residue_map:
                if residue.chain.index != chain_number:
                    continue  # Only modify specified chain
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

        Notes
        -----
        You must have first called findMissingResidues().

        Examples
        --------
        
        Find missing heavy atoms in Abl kinase structure.
        
        >>> fixer = PDBFixer(pdbid='2F4J')
        >>> fixer.findMissingResidues()
        >>> fixer.findMissingAtoms()
        >>> # Retrieve missing atoms.
        >>> missingAtoms = fixer.missingAtoms
        >>> # Retrieve missing terminal atoms.
        >>> missingTerminals = fixer.missingTerminals
        
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
        """Add all missing heavy atoms, as specified by the missingAtoms, missingTerminals, and missingResidues fields.

        Notes
        -----
        You must already have called findMissingAtoms() to have identified atoms to be added.
        
        Examples
        --------
        
        Find missing heavy atoms in Abl kinase structure.
        
        >>> fixer = PDBFixer(pdbid='2F4J')
        >>> fixer.findMissingResidues()
        >>> fixer.findMissingAtoms()
        >>> fixer.addMissingAtoms()

        """
        
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
            if nearest < 0.13:
                
                # Some atoms are very close together.  Run some dynamics while slowly increasing the strength of the
                # repulsive interaction to try to improve the result.
                
                for i in range(10):
                    context.setParameter('C', 0.15*(i+1))
                    integrator.step(200)
                    d = self._findNearestDistance(context, newTopology, newAtoms)
                    if d > nearest:
                        nearest = d
                        state = context.getState(getPositions=True)
                        if nearest >= 0.13:
                            break
                context.setState(state)
                context.setParameter('C', 1.0)
                mm.LocalEnergyMinimizer.minimize(context)
                state = context.getState(getPositions=True)
            
            # Now create a new Topology, including all atoms from the original one and adding the missing atoms.
            
            (newTopology2, newPositions2, newAtoms2, existingAtomMap2) = self._addAtomsToTopology(False, False)
            
            # Copy over the minimized positions for the new atoms.
            
            for a1, a2 in zip(newAtoms, newAtoms2):
                newPositions2[a2.index] = state.getPositions()[a1.index]
            self.topology = newTopology2
            self.positions = newPositions2
    
    def removeHeterogens(self, keepWater=True):
        """Remove all heterogens from the structure.
        
        Parameters
        ----------
        keepWater : bool, optional, default=True
            If True, water molecules will not be removed.

        Examples
        --------

        Remove heterogens in Abl structure complexed with imatinib.

        >>> fixer = PDBFixer(pdbid='2F4J')
        >>> fixer.removeHeterogens(keepWater=False)

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
    
    def addMissingHydrogens(self, pH=7.0):
        """Add missing hydrogen atoms to the structure.
        
        Parameters
        ----------
        pH : float, optional, default=7.0
            The pH based on which to select hydrogens.

        Notes
        -----
        No extensive electrostatic analysis is performed; only default residue pKas are used.

        Examples
        --------

        Examples
        --------

        Add missing hydrogens appropriate for pH 8.

        >>> fixer = PDBFixer(pdbid='1VII')
        >>> fixer.addMissingHydrogens(pH=8.0)

        """
        modeller = app.Modeller(self.topology, self.positions)
        modeller.addHydrogens(pH=pH)
        self.topology = modeller.topology
        self.positions = modeller.positions
    
    def addSolvent(self, boxSize=None, padding=None, positiveIon='Na+', negativeIon='Cl-', ionicStrength=0*unit.molar):
        """Add a solvent box surrounding the structure.
        
        Parameters
        ----------
        boxSize : simtk.openmm.Vec3, optional, default=None
            The size of the box to fill with water.  If specified, padding must not be specified.
        padding : simtk.unit.Quantity compatible with nanometers, optional, default=None
            Padding around macromolecule for filling box with water.  If specified, boxSize must not be specified.
        positiveIon : str, optional, default='Na+'
            The type of positive ion to add.  Allowed values are 'Cs+', 'K+', 'Li+', 'Na+', and 'Rb+'.
        negativeIon : str, optional, default='Cl-'
            The type of negative ion to add.  Allowed values are 'Cl-', 'Br-', 'F-', and 'I-'.
        ionicStrength : simtk.unit.Quantity with units compatible with molar, optional, default=0*molar 
            The total concentration of ions (both positive and negative) to add.  This does not include ions that are added to neutralize the system.
            
        Examples
        --------

        Add missing residues, heavy atoms, and hydrogens, and then solvate with 10 A padding.

        >>> fixer = PDBFixer(pdbid='1VII')
        >>> fixer.findMissingResidues()
        >>> fixer.findMissingAtoms()
        >>> fixer.addMissingAtoms()
        >>> fixer.addMissingHydrogens(pH=8.0)
        >>> fixer.addSolvent(padding=10*unit.angstrom, ionicStrength=0.050*unit.molar)

        """

        modeller = app.Modeller(self.topology, self.positions)
        forcefield = self._createForceField(self.topology, True)
        modeller.addSolvent(forcefield, padding=padding, boxSize=boxSize, positiveIon=positiveIon, negativeIon=negativeIon, ionicStrength=ionicStrength)
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


def main():
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
        fixer = PDBFixer(filename=sys.argv[1])
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
            fixer.addSolvent(boxSize=options.box*unit.nanometer, positiveIon=options.positiveIon, 
                negativeIon=options.negativeIon, ionicStrength=options.ionic*unit.molar)
        with open(options.output, 'w') as f:
            if fixer.source is not None:
                f.write("REMARK   1 PDBFIXER FROM: %s\n" % fixer.source)
            app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

if __name__ == '__main__':
    main()
