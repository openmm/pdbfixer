import pdbfixer
import openmm.app as app
from pathlib import Path

def test_nonstandard():
    """Test adding hydrogens to nonstandard residues."""
    fixer = pdbfixer.PDBFixer(filename=(Path(__file__).parent / "data" / "4JSV.pdb").as_posix())
    fixer.removeChains(chainIndices=[0, 1, 2])
    fixer.addMissingHydrogens()
    for residue in fixer.topology.residues():
        count = sum(1 for atom in residue.atoms() if atom.element.symbol == 'H')
        if residue.name == 'ADP':
            assert count == 15
        if residue.name in ('MG', 'MGF'):
            assert count == 0

def test_leaving_atoms():
    """Test adding hydrogens to a nonstandard residue with leaving atoms."""
    fixer = pdbfixer.PDBFixer(filename=(Path(__file__).parent / "data" / "1BHL.pdb").as_posix())
    fixer.addMissingHydrogens()
    for residue in fixer.topology.residues():
        count = sum(1 for atom in residue.atoms() if atom.element.symbol == 'H')
        if residue.name == 'CAS':
            assert count == 10

def test_registered_template():
    """Test adding hydrogens based on a template registered by the user."""
    fixer = pdbfixer.PDBFixer(filename=(Path(__file__).parent / "data" / "1BHL.pdb").as_posix())

    # Register a template for CAS from which a single hydrogen has been removed.

    pdb = app.PDBFile((Path(__file__).parent / "data" / "CAS.pdb").as_posix())
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.delete([list(modeller.topology.atoms())[-1]])
    terminal = [atom.name in ('H2', 'OXT', 'HXT') for atom in modeller.topology.atoms()]
    fixer.registerTemplate(modeller.topology, modeller.positions, terminal)

    # See if the correct hydrogens get added.

    fixer.addMissingHydrogens()
    for residue in fixer.topology.residues():
        count = sum(1 for atom in residue.atoms() if atom.element.symbol == 'H')
        if residue.name == 'CAS':
            assert count == 9

def test_end_caps():
    """Test adding hydrogens to a chain capped with ACE and NME."""
    fixer = pdbfixer.PDBFixer(filename=(Path(__file__).parent / "data" / "alanine-dipeptide.pdb").as_posix())
    fixer.addMissingHydrogens()
    forcefield = app.ForceField('amber14/protein.ff14SB.xml')
    forcefield.createSystem(fixer.topology)