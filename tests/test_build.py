from Bio.PDB import PDBList
from simtk import unit
from pdbfixer.pdbfixer import *
import os
import sys

# These are tough PDB codes from http://www.umass.edu/microbio/chime/pe_beta/pe/protexpl/badpdbs.htm
pdbcodes = ['1AS5', '1CBN', '1DPO', '1IGY', '1HAG', '1IAO', '4CPA', '1QCQ']

# Set up PDB retrieval.
pdblist = PDBList()

for pdbcode in pdbcodes:
    try:
        input_pdb_filename = pdblist.retrieve_pdb_file(pdbcode, pdir='.')
    except:
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
        sys.exit(1)


# Signal success.
sys.exit(0)

    
