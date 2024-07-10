import pdbfixer
from pathlib import Path
from io import StringIO

def test_nonstandard():
    """Test adding hydrogens to nonstandard residues."""
    content = (Path(__file__).parent / "data" / "4JSV.pdb").read_text()
    fixer = pdbfixer.PDBFixer(pdbfile=StringIO(content))
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
    content = (Path(__file__).parent / "data" / "1BHL.pdb").read_text()
    fixer = pdbfixer.PDBFixer(pdbfile=StringIO(content))
    fixer.addMissingHydrogens()
    for residue in fixer.topology.residues():
        count = sum(1 for atom in residue.atoms() if atom.element.symbol == 'H')
        if residue.name == 'CAS':
            assert count == 10
