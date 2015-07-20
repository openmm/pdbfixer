from nose.tools import ok_, eq_, raises, assert_list_equal
import simtk.openmm.app as app
import simtk.openmm.app.internal.pdbstructure
import pdbfixer
import tempfile
import time

try:
    from urllib.request import urlopen
    from io import StringIO
except:
    from urllib2 import urlopen
    from cStringIO import StringIO

def download_file():
    file = urlopen('http://www.rcsb.org/pdb/files/4JSV.pdb')
    return file.read().decode('utf-8')

def remove_chain_ids_and_verify(file_content, chain_ids_to_remove, expected_chain_ids_remaining):
    # Create a PDBFixer instance for the given pdbid
    fixer = pdbfixer.PDBFixer(pdbfile=StringIO(file_content))
    # Remove specified chains.
    fixer.removeChains(chainIds=chain_ids_to_remove)
    # Check to make sure asserted chains remain.
    chain_ids_remaining = [c.id for c in fixer.topology.chains()]
    assert_list_equal(chain_ids_remaining, expected_chain_ids_remaining)

def test_removechain_ids():
    content = download_file()
    remove_chain_ids_and_verify(content, [], ['B', 'D', 'A', 'C', 'B', 'A'])
    remove_chain_ids_and_verify(content, ['B', 'D'], ['A', 'C', 'A'])
    remove_chain_ids_and_verify(content, ['A', 'C'], ['B', 'D', 'B'])
    remove_chain_ids_and_verify(content, ['B', 'A'], ['D', 'C'])
    remove_chain_ids_and_verify(content, ['B', 'D', 'A', 'C'], [])

def remove_chain_indices_and_verify(file_content, chain_indices_to_remove, expected_chain_ids_remaining):
    # Create a PDBFixer instance for the given pdbid
    fixer = pdbfixer.PDBFixer(pdbfile=StringIO(file_content))
    # Remove specified chains.
    fixer.removeChains(chainIndices=chain_indices_to_remove)
    # Check to make sure asserted chains remain.
    chain_ids_remaining = [c.id for c in fixer.topology.chains()]
    assert_list_equal(chain_ids_remaining, expected_chain_ids_remaining)

def test_removechain_indices():
    content = download_file()
    remove_chain_indices_and_verify(content, [], ['B', 'D', 'A', 'C', 'B', 'A'])
    remove_chain_indices_and_verify(content, [0, 1], ['A', 'C', 'B', 'A'])
    remove_chain_indices_and_verify(content, [2, 3], ['B', 'D', 'B', 'A'])
    remove_chain_indices_and_verify(content, [0, 2], ['D', 'C', 'B', 'A'])
    remove_chain_indices_and_verify(content, [0, 1, 2, 3, 4, 5], [])
