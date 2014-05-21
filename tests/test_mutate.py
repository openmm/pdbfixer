from nose.tools import ok_, eq_, raises
import simtk.openmm.app as app
import pdbfixer
import tempfile

def test_mutate_1():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-GLY"], "A")
    fixer.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, tempfile.TemporaryFile())


def test_mutate_2():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-GLY", "SER-56-ALA"], "A")
    fixer.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, tempfile.TemporaryFile())

@raises(ValueError)
def test_mutate_3_fails():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-GLY", "SER-57-ALA"], "A")
    fixer.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, tempfile.TemporaryFile())


@raises(KeyError)
def test_mutate_4_fails():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-57-WTF", "SER-56-ALA"], "A")
    fixer.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, tempfile.TemporaryFile())


@raises(KeyError)
def test_mutate_5_fails():
    fixer = pdbfixer.PDBFixer(pdbid='1VII')
    fixer.applyMutations(["ALA-1000-GLY", "SER-56-ALA"], "A")
    fixer.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, tempfile.TemporaryFile())
