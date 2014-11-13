from nose.tools import ok_, eq_, raises, assert_items_equal
import simtk.openmm.app as app
import pdbfixer
import tempfile
import copy

def remove_residues_and_verify(pdbid, residueStart, residueEnd, chain_id):
    # Create a PDBFixer instance for the given pdbid
    fixer = pdbfixer.PDBFixer(pdbid=pdbid)
    # Check to make sure asserted chains are removed.
    chains = { c.chain_id : c for c in fixer.structureChains }
    residues_before = copy.deepcopy( chains[chain_id].residues )
    # Remove specified residues.
    fixer.removeResiduesFromChains(residueStart, residueEnd, chainIds=[chain_id])
    # Check to make sure asserted chains are removed.
    chains = { c.chain_id : c for c in fixer.structureChains }
    residues_after = copy.deepcopy( chains[chain_id].residues )

    for residue in residues_after:
        if residue.number in range(residueStart, residueEnd+1):
            raise Exception("Found residue %d, but should have removed from %d to %d in chain %s" % (residue.number, residueStart, residueEnd, chain_id))

def test_removeresidues():
    remove_residues_and_verify('1VII', 53, 67, 'A')
    remove_residues_and_verify('4J6F', 1385, 2001, 'B')
