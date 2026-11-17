import io
import pdbfixer

# SEQRES starts two residues before the first resolved residue (ALA 1), so the
# residues built at the N terminus should be numbered -1 and 0.
PDB = """\
SEQRES   1 A    5  MET GLY ALA VAL GLU
ATOM      1  N   ALA A   1      -1.200   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       1.200   0.400   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.400   1.600   0.000  1.00 20.00           O
ATOM      5  N   VAL A   2       2.600   0.000   0.000  1.00 20.00           N
ATOM      6  CA  VAL A   2       3.800   0.000   0.000  1.00 20.00           C
ATOM      7  C   VAL A   2       5.000   0.400   0.000  1.00 20.00           C
ATOM      8  O   VAL A   2       5.200   1.600   0.000  1.00 20.00           O
ATOM      9  N   GLU A   3       6.400   0.000   0.000  1.00 20.00           N
ATOM     10  CA  GLU A   3       7.600   0.000   0.000  1.00 20.00           C
ATOM     11  C   GLU A   3       8.800   0.400   0.000  1.00 20.00           C
ATOM     12  O   GLU A   3       9.000   1.600   0.000  1.00 20.00           O
TER      13      GLU A   3
END
"""

def test_negative_residue_ids():
    fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(PDB))
    fixer.findMissingResidues()
    assert fixer.missingResidues == {(0, 0): ['MET', 'GLY']}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    ids = [residue.id for residue in fixer.topology.residues()]
    assert ids == ['-1', '0', '1', '2', '3']
