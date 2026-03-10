from pathlib import Path

from openmm import app, unit

import pdbfixer


def test_findOXTPosition():
    """Test that OXT is added at the correct position."""
    fixer = pdbfixer.PDBFixer(filename=(Path(__file__).parent / "data" / "1BHL.pdb").as_posix())

    # Record the original position of OXT, then delete it.

    originalPos = [pos for pos, atom in zip(fixer.positions, fixer.topology.atoms()) if atom.name == "OXT"]
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.delete([atom for atom in modeller.topology.atoms() if atom.name == "OXT"])
    fixer.topology = modeller.topology
    fixer.positions = modeller.positions

    # Have PDBFixer add it back and make sure it is sufficiently close to the original position.

    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    newPos = [pos for pos, atom in zip(fixer.positions, fixer.topology.atoms()) if atom.name == "OXT"]
    assert unit.norm(newPos[0] - originalPos[0]) < 0.01 * unit.nanometer
