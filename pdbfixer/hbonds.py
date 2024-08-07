import openmm as mm
import openmm.app as app
import openmm.unit as unit
from openmm.app.element import hydrogen, oxygen, nitrogen, fluorine
import numpy as np
import itertools
from collections import defaultdict
from collections.abc import Iterator
from typing import List

class HydrogenBondOptimizer(object):
    
    def __init__(self, fixer):
        topology, positions = self.protonateHistidines(fixer)

        # Record the parent atom that every hydrogen is bonded to.

        parent = {}
        for atom1, atom2 in topology.bonds():
            if atom1.element == hydrogen and atom2.element != hydrogen:
                parent[atom1] = atom2
            elif atom2.element == hydrogen and atom1.element != hydrogen:
                parent[atom2] = atom1

        # Record the donors and acceptors in each residue.

        self.residueDonors = {}
        self.residueAcceptors = {}
        for residue in topology.residues():
            self.residueDonors[residue] = []
            self.residueAcceptors[residue] = []
            for atom in residue.atoms():
                if atom.element == hydrogen and parent[atom].element in (oxygen, nitrogen, fluorine):
                    self.residueDonors[residue].append((atom, parent[atom]))
                elif atom.element in (oxygen, nitrogen, fluorine):
                    self.residueAcceptors[residue].append(atom)

        # Record the groups of atoms we will try to rotate.

        rotations = {
            'ASN': (('CB', 'CG'), ('OD1', 'ND2', 'HD21', 'HD22')),
            'GLN': (('CG', 'CD'), ('OE1', 'NE2', 'HE21', 'HE22')),
            'HIS': (('CB', 'CG'), ('ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HD2', 'HE1', 'HE2'))
        }
        residueRotations = {}
        for residue in topology.residues():
            if residue.name in rotations:
                axis, atoms = rotations[residue.name]
                residueRotations[residue] = ([atom for atom in residue.atoms() if atom.name in axis],
                                             [atom for atom in residue.atoms() if atom.name in atoms])

        # For each residue we will be modifying, record all other residues that are close enough
        # to potentially form hydrogen bonds.

        self.residueNeighbors = {}
        residueBounds = dict((residue, BoundingSphere(residue.atoms(), positions)) for residue in topology.residues())
        for res1 in residueRotations:
            self.residueNeighbors[res1] = []
            bounds = BoundingSphere(residueRotations[res1][1], positions)
            for res2 in topology.residues():
                if res1 != res2 and bounds.distance(residueBounds[res2]) < 0.5:
                    close = False
                    for atom1 in residueRotations[res1][1]:
                        for atom2 in res2.atoms():
                            delta = positions[atom1.index]-positions[atom2.index]
                            if unit.norm(delta) < 0.5:
                                close = True
                                break
                        if close:
                            break
                    if close:
                        self.residueNeighbors[res1].append(res2)

        # Divide them into clusters that need to be analyzed together.

        clusters = [[res] for res in self.residueNeighbors]
        converged = False
        while not converged:
            converged = True
            for i in range(len(clusters)):
                if i >= len(clusters):
                    break
                c1 = clusters[i]
                for j in range(i+1, len(clusters)):
                    c2 = clusters[j]
                    if any(res2 in self.residueNeighbors[res1] for res1, res2 in itertools.product(c1, c2)):
                        c1 += c2
                        del clusters[j]
                        converged = False
                        break

        # Create the rotated version of the coordinates.

        rotatedPositions = positions[:]
        for residue in residueRotations:
            axis, atoms = residueRotations[residue]
            center = positions[axis[1].index]
            e = center-positions[axis[0].index]
            e /= unit.norm(e)
            for atom in atoms:
                d = positions[atom.index]-center
                dot = np.dot(e, d)
                rotatedPositions[atom.index] = center - d + 2*dot*e

        # For each variable residue, consider hydrogen bonds to other residues that are
        # *not* variable.  How does the number change?

        residueRotationHbondChange = defaultdict(int)
        for residue in self.residueNeighbors:
            for neighbor in self.residueNeighbors[residue]:
                if neighbor not in self.residueNeighbors:
                    before = self.countHbonds(residue, neighbor, positions, positions)
                    after = self.countHbonds(residue, neighbor, rotatedPositions, positions)
                    residueRotationHbondChange[residue] += after-before

        # Now consider pairs of residues that are close enough to form hydrogen bonds to
        # each other.

        pairRotationHbonds = {}
        for residue in self.residueNeighbors:
            for neighbor in self.residueNeighbors[residue]:
                if neighbor.index > residue.index and neighbor in self.residueNeighbors:
                    pairRotationHbonds[(residue, neighbor)] = [
                        [
                            self.countHbonds(residue, neighbor, positions, positions),
                            self.countHbonds(residue, neighbor, positions, rotatedPositions)
                        ],
                        [
                            self.countHbonds(residue, neighbor, rotatedPositions, positions),
                            self.countHbonds(residue, neighbor, rotatedPositions, rotatedPositions)
                        ]
                    ]

        # Optimize each cluster.

        newPositions = positions[:]
        for cluster in clusters:
            bestCombination = None
            bestBonds = -1

            # Loop over all possible combinations of which residues to rotate.

            for combination in itertools.product(*[[0,1]]*len(cluster)):
                bonds = 0
                for rotated, residue in zip(combination, cluster):
                    if rotated == 1:
                        bonds += residueRotationHbondChange[residue]
                    for neighbor in self.residueNeighbors[residue]:
                        if (residue, neighbor) in pairRotationHbonds:
                            neighborRotated = combination[cluster.index(neighbor)]
                            bonds += pairRotationHbonds[(residue, neighbor)][rotated][neighborRotated]
                if bonds > bestBonds:
                    bestCombination = combination
                    bestBonds = bonds

            # Set the positions to the optimal set of rotations.

            for rotated, residue in zip(bestCombination, cluster):
                if rotated == 1:
                    for atom in residueRotations[residue][1]:
                        newPositions[atom.index] = rotatedPositions[atom.index]

        # Create the final topology and positions.

        fixer.topology, fixer.positions = self.deleteExtraHydrogens(fixer, topology, positions)

    def protonateHistidines(self, fixer):
        """Add a second hydrogen to neutral histidines so we can identify bonds formed by both of them."""
        positions = fixer.positions.value_in_unit(unit.nanometer)
        newTopology = app.Topology()
        newTopology.setPeriodicBoxVectors(fixer.topology.getPeriodicBoxVectors())
        newPositions = []
        newAtoms = {}
        for chain in fixer.topology.chains():
            newChain = newTopology.addChain(chain.id)
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                for atom in residue.atoms():
                    newAtoms[atom] = newTopology.addAtom(atom.name, atom.element, newResidue)
                    newPositions.append(positions[atom.index])
                if residue.name == 'HIS':
                    atomsByName = dict((atom.name, atom) for atom in newResidue.atoms())
                    posByName = dict((atom.name, newPositions[atom.index]) for atom in newResidue.atoms())
                    if all(name in atomsByName for name in ('CG', 'ND1', 'CD2', 'CE1', 'NE2')):
                        # If this residue has one of the two hydrogens, add the other.
                        if 'HD1' in atomsByName and 'HE2' not in atomsByName:
                            newAtom = newTopology.addAtom('HE2', hydrogen, newResidue)
                            newTopology.addBond(newAtom, atomsByName['NE2'])
                            axis = posByName['NE2'] - 0.5*(posByName['CD2']-posByName['CE1'])
                            axis /= unit.norm(axis)
                            newPositions.append(posByName['NE2'] + axis)
                        elif 'HD1' not in atomsByName and 'HE2' in atomsByName:
                            newAtom = newTopology.addAtom('HD1', hydrogen, newResidue)
                            newTopology.addBond(newAtom, atomsByName['ND1'])
                            axis = posByName['ND1'] - 0.5*(posByName['CG']-posByName['CE1'])
                            axis /= unit.norm(axis)
                            newPositions.append(posByName['ND1'] + axis)
        for atom1, atom2 in fixer.topology.bonds():
            newTopology.addBond(newAtoms[atom1], newAtoms[atom2])
        return newTopology, newPositions

    def deleteExtraHydrogens(self, fixer, topology, positions):
        """We may have added an extra hydrogen to some HIS residues.  Figure out which one to keep
        and delete the other one."""
        toDelete = []
        for oldRes, newRes in zip(fixer.topology.residues(), topology.residues()):
            if oldRes.name == 'HIS' and len(oldRes) != len(newRes):
                # See whether each of the hydrogens forms a hydrogen bond.

                atomsByName = dict((atom.name, atom) for atom in newRes.atoms())
                hd1 = atomsByName['HD1']
                nd1 = atomsByName['ND1']
                he2 = atomsByName['HE2']
                ne2 = atomsByName['NE2']
                hd1Acceptor = None
                he2Acceptor = None
                for neighbor in self.residueNeighbors[newRes]:
                    for a in self.residueAcceptors[neighbor]:
                        if self.isHbond(positions[nd1.index], positions[hd1.index], positions[a.index]):
                            hd1Acceptor = neighbor
                        if self.isHbond(positions[ne2.index], positions[he2.index], positions[a.index]):
                            he2Acceptor = neighbor
                if hd1Acceptor is not None and he2Acceptor is None:
                    # HD1 forms a hydrogen bond.  Delete HE2.
                    toDelete.append(he2)
                elif hd1Acceptor is None and he2Acceptor is not None:
                    # HE2 forms a hydrogen bond.  Delete HD1.
                    toDelete.append(hd1)
                elif hd1Acceptor is None and he2Acceptor is None:
                    # Neither one forms a hydrogen bond.  Delete the one that was not present
                    # in the original Topology.
                    if any(atom.name == 'HD1' for atom in oldRes.atoms()):
                        toDelete.append(he2)
                    else:
                        toDelete.append(hd1)
                else:
                    # Both form hydrogen bonds.  If one is bonded to a rotatable residue, keep that one.
                    # That residue was optimized based on the assumption it could form a hydrogen bond
                    # to this one.
                    if hd1Acceptor.name in ('ASN', 'GLN', 'HIS'):
                        toDelete.append(he2)
                    elif he2Acceptor.name in ('ASN', 'GLN', 'HIS'):
                        toDelete.append(hd1)
                    else:
                        # Delete the one that was not present in the original Topology.
                        if any(atom.name == 'HD1' for atom in oldRes.atoms()):
                            toDelete.append(he2)
                        else:
                            toDelete.append(hd1)

        # Create the new Topology and positions.

        modeller = app.Modeller(topology, positions)
        modeller.delete(toDelete)
        return modeller.topology, modeller.positions

    def isHbond(self, d: mm.Vec3, h: mm.Vec3, a: mm.Vec3):
        """Decide whether three atoms could form a hydrogen bond."""
        if unit.norm(d-a) > 0.35:
            return False
        deltaDH = h-d
        deltaHA = a-h
        deltaDH /= unit.norm(deltaDH)
        deltaHA /= unit.norm(deltaHA)
        return np.arccos(np.dot(deltaDH, deltaHA)) < 50*np.pi/180

    def countHbonds(self, res1: app.Residue, res2: app.Residue, pos1: List[mm.Vec3], pos2: List[mm.Vec3]):
        """Count the number of hydrogen bonds between two residues."""
        count = 0
        for h, d in self.residueDonors[res1]:
            for a in self.residueAcceptors[res2]:
                if self.isHbond(pos1[d.index], pos1[h.index], pos2[a.index]):
                    count += 1
        for h, d in self.residueDonors[res2]:
            for a in self.residueAcceptors[res1]:
                if self.isHbond(pos2[d.index], pos2[h.index], pos1[a.index]):
                    count += 1
        return count


class BoundingSphere(object):
    """Computes a bounding sphere for a set of atoms.  This is used to accelerate searches for neighboring residues."""

    def __init__(self, atoms: Iterator[app.Atom], positions: List[mm.Vec3]):
        atoms = list(atoms)
        minRange = mm.Vec3(*(min((positions[atom.index][i] for atom in atoms)) for i in range(3)))
        maxRange = mm.Vec3(*(max((positions[atom.index][i] for atom in atoms)) for i in range(3)))
        self.center = 0.5*(minRange+maxRange)
        self.radius = max(unit.norm(self.center-positions[atom.index]) for atom in atoms)

    def distance(self, other: 'BoundingSphere'):
        return max(0, unit.norm(self.center-other.center) - self.radius - other.radius)
