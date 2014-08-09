from nose.tools import ok_, eq_, raises
import simtk.openmm.app as app
import pdbfixer
import tempfile

def remove_chainids_and_verify(pdbid, chain_ids_to_remove, expected_chain_ids_remaining):
    # Create a PDBFixer instance for the given pdbid
    fixer = pdbfixer.PDBFixer(pdbid=pdbid)
    # Remove specified chains.
    fixer.removeChains(chainIds=chain_ids_to_remove)
    # Check to make sure asserted chains remain.
    chain_ids_remaining = [ chain.id for chain in fixer.topology.chains() ]
    assertItemsEqual(chain_ids_remaining, expected_chain_ids_remaining)

def test_removechainids():
    remove_chainids_and_verify('4JSV', [], ['B', 'D', 'A', 'C'])
    remove_chainids_and_verify('4JSV', ['B', 'D'], ['A', 'C'])
    remove_chainids_and_verify('4JSV', ['A', 'C'], ['B', 'D'])
    remove_chainids_and_verify('4JSV', ['B', 'A'], ['D', 'C'])
    remove_chainids_and_verify('4JSV', ['B', 'D', 'A', 'C'], [])
