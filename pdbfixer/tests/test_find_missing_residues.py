import os
import warnings
from pathlib import Path

import pytest

import pdbfixer


def createFixer(filename):
    return pdbfixer.PDBFixer(filename=(Path(__file__).parent / "data" / filename).as_posix())


def missingResiduesByChainId(fixer):
    """Group the contents of missingResidues by chain ID, in the order they appear in the sequence."""
    chains = list(fixer.topology.chains())
    result = {}
    for key in sorted(fixer.missingResidues):
        result.setdefault(chains[key[0]].id, []).extend(fixer.missingResidues[key])
    return result


def test_insertionCodesAndSkippedNumbers():
    """Test a chain whose numbering both uses insertion codes and skips numbers.

    The chain has a genuine two residue gap between residues 52 and 55, and its numbering also
    jumps from 57 to 60 without any residue being absent.  Only the former should be reported.
    """
    fixer = createFixer("insertion_codes.pdb")
    fixer.findMissingResidues()
    assert fixer.missingResidues == {(0, 2): ["ASP", "GLU"]}


def test_insertionCodesAreNotTreatedAsGaps():
    """Test that the missing residues of such a chain, and only those, are built."""
    fixer = createFixer("insertion_codes.pdb")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    residues = [r.name for r in fixer.topology.residues()]
    assert residues == ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
                        "MET", "ASN", "PRO", "GLN"]


def test_warnsWhenChainCannotBeAligned():
    """Test that a chain we fail to align produces a warning rather than being silently skipped."""
    fixer = createFixer("unalignable_order.pdb")
    with pytest.warns(UserWarning, match="Unable to align chain A"):
        fixer.findMissingResidues()
    assert fixer.missingResidues == {}


def test_noWarningForHeterogens():
    """Test that heterogens and waters sharing a chain ID with a polymer do not produce warnings."""
    fixer = createFixer("4JSV.pdb")
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        fixer.findMissingResidues()
    assert [str(r.message) for r in record if issubclass(r.category, UserWarning)] == []


@pytest.mark.skipif(os.getenv("GITHUB_ACTION") is not None, reason="Downloads files from RCSB")
@pytest.mark.parametrize(
    "pdbid,expected",
    [
        # An intact IgG whose heavy chains are Kabat numbered.  Nothing is missing from any chain,
        # although the numbering of chains B and D skips 23 numbers and uses 8 insertion codes.
        ("1IGT", {}),
        # An Fab whose heavy chain is Kabat numbered.  These are the residues REMARK 465 lists.
        ("1A0Q", {"L": ["ASP"],
                  "H": ["GLU", "GLY", "ARG", "SER", "ASN",
                        "GLY", "SER", "ALA", "ALA", "GLN", "THR", "ASN"]}),
        # A structure with sequential numbering, included to catch regressions.
        ("1VII", {}),
    ],
)
def test_missingResiduesMatchDeposition(pdbid, expected):
    """Test that the missing residues found are the ones the deposition lists in REMARK 465."""
    fixer = pdbfixer.PDBFixer(pdbid=pdbid)
    fixer.findMissingResidues()
    assert missingResiduesByChainId(fixer) == expected
