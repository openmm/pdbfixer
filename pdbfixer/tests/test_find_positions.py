from pathlib import Path

import numpy as np
import pytest
from openmm import Vec3, app, unit

import pdbfixer
from pdbfixer.pdbfixer import _overlayPoints


def _load1BHL():
    """Load the 1BHL test structure."""
    return pdbfixer.PDBFixer(filename=(Path(__file__).parent / "data" / "1BHL.pdb").as_posix())


def _deleteAndRebuild(fixer, shouldDelete):
    """Delete the atoms selected by shouldDelete, then have PDBFixer add them back without building missing residues."""
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.delete([atom for atom in modeller.topology.atoms() if shouldDelete(atom)])
    fixer.topology = modeller.topology
    fixer.positions = modeller.positions
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()


def _oxtPositions(fixer):
    """Return the positions of all OXT atoms."""
    return [pos for pos, atom in zip(fixer.positions, fixer.topology.atoms()) if atom.name == "OXT"]


def _residueHandedness(fixer, shouldInclude):
    """Return the handedness (+1 or -1) of the selected alpha carbons, in residue order."""
    positions = fixer.positions.value_in_unit(unit.nanometer)
    handedness = []
    for residue in fixer.topology.residues():
        if shouldInclude(residue):
            pos = {atom.name: positions[atom.index] for atom in residue.atoms()}
            toN = pos["N"]-pos["CA"]
            toC = pos["C"]-pos["CA"]
            toCB = pos["CB"]-pos["CA"]
            handedness.append(np.sign(np.dot(toN, np.cross(toC, toCB))))
    return handedness


def test_overlayPointsNeverReflects():
    """Test that _overlayPoints returns a proper rotation even when a reflection would fit better."""
    points2 = [Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(1, 1, 0.1)]
    points1 = [Vec3(p[0], p[1], -p[2]) for p in points2]
    (_, rotate, _) = _overlayPoints(points1, points2)
    assert np.isclose(np.linalg.det(rotate), 1)


def test_findOXTPosition():
    """Test that OXT is added at the correct position."""
    fixer = _load1BHL()

    # Record the original position of OXT, then delete it and have PDBFixer add it back.

    originalPos = _oxtPositions(fixer)
    _deleteAndRebuild(fixer, lambda atom: atom.name == "OXT")

    # Make sure it is sufficiently close to the original position.

    assert unit.norm(_oxtPositions(fixer)[0] - originalPos[0]) < 0.01 * unit.nanometer


def test_addedCBChirality():
    """Test that missing CB atoms are added with the correct chirality."""
    fixer = _load1BHL()

    # Record the handedness of every alanine, then delete its CB and have PDBFixer add it back.

    originalHandedness = _residueHandedness(fixer, lambda residue: residue.name == "ALA")
    assert len(originalHandedness) > 0
    _deleteAndRebuild(fixer, lambda atom: atom.residue.name == "ALA" and atom.name == "CB")

    # Make sure none of them was built as its mirror image.  Adding atoms keeps the residues
    # in the same order, so the lists can be compared directly.

    assert _residueHandedness(fixer, lambda residue: residue.name == "ALA") == originalHandedness


@pytest.mark.parametrize("residueId", ["57", "72", "112"])
def test_addedSideChainChirality(residueId):
    """Test that the backbone oxygen does not put a rebuilt side chain on the wrong side."""
    fixer = _load1BHL()

    def isSelected(residue):
        return residue.chain.id == "A" and residue.id == residueId

    # SER57, VAL72, and THR112 have backbone oxygen positions that tilt the template
    # fit enough to put CB on the wrong side, even when the fit is a proper rotation.

    originalHandedness = _residueHandedness(fixer, isSelected)
    assert len(originalHandedness) == 1
    _deleteAndRebuild(fixer, lambda atom: isSelected(atom.residue) and atom.name not in ("N", "CA", "C", "O", "OXT"))

    assert _residueHandedness(fixer, isSelected) == originalHandedness
