import openmm.app as app
import pdbfixer
from pytest import raises
from pathlib import Path

def test_mutate_1():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-GLY"], "A")
    fixer.findMissingResidues()     
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()        
    fixer.addMissingHydrogens(7.0)

    new_residue57 = list(fixer.topology.residues())[16]
    assert new_residue57.name == "GLY", "Name of mutated residue did not change correctly!"
    assert len(list(new_residue57.atoms())) == 7, "Should have 7 atoms in GLY 56"
    
    atom_names = set([atom.name for atom in new_residue57.atoms()])
    desired_atom_names = set(["N", "CA", "C", "O", "H", "HA3", "HA2"])
    assert atom_names == desired_atom_names, "Atom Names did not match for GLY 56"    


def test_mutate_2():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-LEU", "SER-56-ALA"], "A")
    fixer.findMissingResidues()     
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()        
    fixer.addMissingHydrogens(7.0)
    new_residue57 = list(fixer.topology.residues())[16]
    new_residue56 = list(fixer.topology.residues())[15]
    assert new_residue57.name == "LEU", "Name of mutated residue did not change correctly!"
    assert new_residue56.name == "ALA", "Name of mutated residue did not change correctly!"
    
    assert len(list(new_residue56.atoms())) == 10, "Should have 10 atoms in ALA 56"
    assert len(list(new_residue57.atoms())) == 19, "Should have 19 atoms in LEU 57"
    
    atom_names = set([atom.name for atom in new_residue56.atoms()])
    desired_atom_names = set(["N", "CA", "CB", "C", "O", "H", "HA", "HB1", "HB2", "HB3"])
    assert atom_names == desired_atom_names, "Atom Names did not match for ALA 56"

    atom_names = set([atom.name for atom in new_residue57.atoms()])
    desired_atom_names = set(["C", "N", "CA", "CB", "CG", "CD1", "CD2", "O", "H", "HA", "HB2", "HB3", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23", "HG"])
    assert atom_names == desired_atom_names, "Atom Names did not match for LEU 57"

def test_mutate_3_fails():
    with raises(ValueError):
        fixer = pdbfixer.PDBFixer(pdbid='1VII')
        fixer.applyMutations(["ALA-57-GLY", "SER-57-ALA"], "A")

def test_mutate_4_fails():
    with raises(KeyError):
        fixer = pdbfixer.PDBFixer(pdbid='1VII')
        fixer.applyMutations(["ALA-1000-GLY", "SER-56-ALA"], "A")

def test_mutate_multiple_copies_of_chain_A():
    fixer = pdbfixer.PDBFixer(pdbid='1OL5')
    fixer.applyMutations(['TPO-287-THR','TPO-288-THR'], "A")

def test_mutate_to_nonstandard():
    """Test mutating to a nonstandard residue defined with registerTemplate()."""
    fixer = pdbfixer.PDBFixer(filename=(Path(__file__).parent / "data" / "1BHL.pdb").as_posix())
    pdb = app.PDBFile((Path(__file__).parent / "data" / "CAS.pdb").as_posix())
    terminal = [atom.name in ('H2', 'OXT', 'HXT') for atom in pdb.topology.atoms()]
    fixer.registerTemplate(pdb.topology, pdb.positions, terminal)
    fixer.applyMutations(["SER-57-CAS", "ILE-60-CAS", "ASP-207-CAS"], "A")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    residues = list(fixer.topology.residues())
    for i in (0, 3, 134):
        assert residues[i].name == "CAS"
        atoms = list(residues[i].atoms())
        assert sum(1 for a in atoms if a.name == 'H2') == (1 if i == 0 else 0)
        assert sum(1 for a in atoms if a.name == 'OXT') == (1 if i == 134 else 0)
        assert sum(1 for a in atoms if a.name == 'HXT') == (1 if i == 134 else 0)

def test_download_template():
    """Test mutating to a nonstandard residue defined in the PDB."""
    fixer = pdbfixer.PDBFixer(filename=(Path(__file__).parent / "data" / "1BHL.pdb").as_posix())
    fixer.applyMutations(["SER-57-SEP", "ILE-60-SEP", "ASP-207-SEP"], "A")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    residues = list(fixer.topology.residues())
    for i in (0, 3, 134):
        assert residues[i].name == "SEP"
        atoms = list(residues[i].atoms())
        assert sum(1 for a in atoms if a.element == app.element.phosphorus) == 1
        assert sum(1 for a in atoms if a.name == 'OXT') == (1 if i == 134 else 0)

    # Check a few bonds to make sure the mutated residue has the ones it's supposed to.

    bonds = list(residues[3].bonds())
    assert sum(1 for a1, a2 in bonds if {a1.name, a2.name} == {'N', 'CA'}) == 1
    assert sum(1 for a1, a2 in bonds if {a1.name, a2.name} == {'CB', 'OG'}) == 1
    assert sum(1 for a1, a2 in bonds if {a1.name, a2.name} == {'P', 'OG'}) == 1
    assert sum(1 for a1, a2 in bonds if {a1.name, a2.name} == {'P', 'O2P'}) == 1