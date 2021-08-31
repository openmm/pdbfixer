import openmm.app as app
import pdbfixer
import tempfile
from pytest import raises

def test_mutate_1():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-GLY"], "A")
    fixer.findMissingResidues()     
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()        
    fixer.addMissingHydrogens(7.0)
    with tempfile.NamedTemporaryFile(mode='w+') as temp_pdb:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, temp_pdb)
        temp_pdb.flush()
        pdb = app.PDBFile(temp_pdb.name)
    
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
    temp_pdb = tempfile.NamedTemporaryFile(mode='w+')
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
        fixer.applyMutations(["ALA-57-WTF", "SER-56-ALA"], "A")


def test_mutate_5_fails():
    with raises(KeyError):
        fixer = pdbfixer.PDBFixer(pdbid='1VII')
        fixer.applyMutations(["ALA-1000-GLY", "SER-56-ALA"], "A")

def test_mutate_multiple_copies_of_chain_A():
    fixer = pdbfixer.PDBFixer(pdbid='1OL5')
    fixer.applyMutations(['TPO-287-THR','TPO-288-THR'], "A")

