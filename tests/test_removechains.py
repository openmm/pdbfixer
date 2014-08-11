from nose.tools import ok_, eq_, raises, assert_items_equal
import simtk.openmm.app as app
import pdbfixer
import tempfile

def remove_chain_ids_and_verify(pdbid, chain_ids_to_remove, expected_chain_ids_remaining):
    # Create a PDBFixer instance for the given pdbid
    fixer = pdbfixer.PDBFixer(pdbid=pdbid)
    # Remove specified chains.
    fixer.removeChains(chainIds=chain_ids_to_remove)
    # Check to make sure asserted chains remain.
    chain_ids_remaining = [c.chain_id for c in fixer.structureChains]
    assert_items_equal(chain_ids_remaining, expected_chain_ids_remaining)

def test_removechain_ids():
    remove_chain_ids_and_verify('4JSV', [], ['B', 'D', 'A', 'C', 'B', 'A'])
    remove_chain_ids_and_verify('4JSV', ['B', 'D'], ['A', 'C', 'A'])
    remove_chain_ids_and_verify('4JSV', ['A', 'C'], ['B', 'D', 'B'])
    remove_chain_ids_and_verify('4JSV', ['B', 'A'], ['D', 'C'])
    remove_chain_ids_and_verify('4JSV', ['B', 'D', 'A', 'C'], [])

def remove_chain_indices_and_verify(pdbid, chain_indices_to_remove, expected_chain_ids_remaining):
    # Create a PDBFixer instance for the given pdbid
    fixer = pdbfixer.PDBFixer(pdbid=pdbid)
    # Remove specified chains.
    fixer.removeChains(chainIndices=chain_indices_to_remove)
    # Check to make sure asserted chains remain.
    chain_ids_remaining = [c.chain_id for c in fixer.structureChains]
    assert_items_equal(chain_ids_remaining, expected_chain_ids_remaining)

def test_removechain_indices():
    remove_chain_indices_and_verify('4JSV', [], ['B', 'D', 'A', 'C', 'B', 'A'])
    remove_chain_indices_and_verify('4JSV', [0, 1], ['A', 'C', 'B', 'A'])
    remove_chain_indices_and_verify('4JSV', [2, 3], ['B', 'D', 'B', 'A'])
    remove_chain_indices_and_verify('4JSV', [0, 2], ['D', 'C', 'B', 'A'])
    remove_chain_indices_and_verify('4JSV', [0, 1, 2, 3, 4, 5], [])
