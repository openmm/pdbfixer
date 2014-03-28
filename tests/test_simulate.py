from Bio.PDB import PDBList
from simtk import unit
from pdbfixer.pdbfixer import *
import os
import sys

def test_build():
    # These are tough PDB codes from http://www.umass.edu/microbio/chime/pe_beta/pe/protexpl/badpdbs.htm
    pdbcodes = ['1AS5', '1CBN', '1DPO', '1IGY', '1HAG', '1IAO', '4CPA', '1QCQ']
    pdbcodes = ['1VII'] # DEBUG: just use one

    # Set up PDB retrieval.
    pdblist = PDBList(server=PDBList.alternative_download_url)

    success = True

    for pdbcode in pdbcodes:
        print pdbcode

        try:
            print "Attempting to retrieve PDB code '%s' from %s..." % (pdbcode, PDBList.alternative_download_url)
            input_pdb_filename = pdblist.retrieve_pdb_file(pdbcode, pdir='.')
        except Exception as e:
            print str(e)
            print "Could not download PDB code '%s'" % pdbcode
            continue

        output_pdb_filename = 'output.pdb'

        # PDB setup parameters.
        # TODO: Try several combinations?
        pH = 7.0
        ionic = 50.0 * unit.millimolar
        box = 10.0 * unit.angstrom
        positiveIon = 'Na+'
        negativeIon = 'Cl-'

        infile = open(input_pdb_filename)
        outfile = open(output_pdb_filename, 'w')

        try:
            fixer = PDBFixer(PdbStructure(infile))
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.removeHeterogens(False)
            fixer.addMissingHydrogens(pH)
            #fixer.addSolvent(box*unit.nanometer, positiveIon, negativeIon, ionic*unit.molar)
            app.PDBFile.writeFile(fixer.topology, fixer.positions, outfile)
            infile.close()
            outfile.close()
            
            # Delete input file.
            os.remove(input_pdb_filename)
            os.remove(output_pdb_filename)

        except Exception as e:
            print str(e)
            success = False

    if not success:
        raise Exception("build test failed on one or more PDB files.")

