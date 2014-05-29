from nose.tools import ok_, eq_, raises
import simtk.openmm.app as app
import pdbfixer
import tempfile

def test_mutate_1():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-GLY"], "A")
    fixer.findMissingResidues()     
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()        
    fixer.addMissingHydrogens(7.0)
    temp_pdb = tempfile.NamedTemporaryFile()
    app.PDBFile.writeFile(fixer.topology, fixer.positions, temp_pdb)
    pdb = app.PDBFile(temp_pdb.name)
    
    new_residue57 = list(fixer.topology.residues())[16]
    assert new_residue57.name == "GLY", "Name of mutated residue did not change correctly!"
    assert len(list(new_residue57.atoms())) == 7, "Should have 7 atoms in GLY 56"    


def test_mutate_2():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-LEU", "SER-56-ALA"], "A")
    fixer.findMissingResidues()     
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()        
    fixer.addMissingHydrogens(7.0)
    temp_pdb = tempfile.NamedTemporaryFile()
    new_residue57 = list(fixer.topology.residues())[16]
    new_residue56 = list(fixer.topology.residues())[15]
    assert new_residue57.name == "LEU", "Name of mutated residue did not change correctly!"
    assert new_residue56.name == "ALA", "Name of mutated residue did not change correctly!"
    
    assert len(list(new_residue56.atoms())) == 10, "Should have 10 atoms in ALA 56"
    assert len(list(new_residue57.atoms())) == 19, "Should have 19 atoms in LEU 57"


@raises(ValueError)
def test_mutate_3_fails():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-GLY", "SER-57-ALA"], "A")

@raises(KeyError)
def test_mutate_4_fails():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-WTF", "SER-56-ALA"], "A")


@raises(KeyError)
def test_mutate_5_fails():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-1000-GLY", "SER-56-ALA"], "A")

