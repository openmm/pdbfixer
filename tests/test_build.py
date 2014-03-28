from Bio.PDB import PDBList
from simtk import unit
from pdbfixer.pdbfixer import *
import os
import sys

def test_build():
    # These are tough PDB codes from http://www.umass.edu/microbio/chime/pe_beta/pe/protexpl/badpdbs.htm
    pdbcodes = ['1AS5', '1CBN', '1DPO', '1IGY', '1HAG', '1IAO', '4CPA', '1QCQ']
    pdbcodes = ['1VII'] # DEBUG

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

        print "Running PDBFixer..."
        infile = open(input_pdb_filename)
        outfile = open(output_pdb_filename, 'w')

        try:
            print "Creating PDBFixer..."
            fixer = PDBFixer(PdbStructure(infile))
            print "Finding missing residues..."
            fixer.findMissingResidues()
            print "Finding nonstandard residues..."
            fixer.findNonstandardResidues()
            print "Replacing nonstandard residues..."
            fixer.replaceNonstandardResidues()
            print "Finding missing atoms..."
            fixer.findMissingAtoms()
            print "Adding missing atoms..."
            fixer.addMissingAtoms()
            print "Removing heterogens..."
            fixer.removeHeterogens(False)
            print "Adding missing hydrogens..."
            fixer.addMissingHydrogens(pH)
            #fixer.addSolvent(box*unit.nanometer, positiveIon, negativeIon, ionic*unit.molar)
            print "Writing PDB file..."
            app.PDBFile.writeFile(fixer.topology, fixer.positions, outfile)
            infile.close()
            outfile.close()
            print "Done."
            
            # Delete input file.
            os.remove(input_pdb_filename)
            os.remove(output_pdb_filename)

        except Exception as e:
            print str(e)
            success = False

    if not success:
        raise Exception("build test failed on one or more PDB files.")
            
test_build()
