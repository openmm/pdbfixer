
import pdbfixer
import simtk.openmm 
import Bio.PDB

import os
import sys
import numpy

def simulate(pdbcode, pdb_filename):
    from simtk.openmm import app
    import simtk.openmm as mm
    from simtk import unit
    from sys import stdout

    # Load the PDB file.
    pdb = app.PDBFile(pdb_filename)
    
    # Set up implicit solvent forcefield.
    forcefield = app.ForceField('amber99sbildn.xml')
    
    # Create the system.
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)

    # Create an integrator.
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 91.0/unit.picoseconds, 1.0*unit.femtoseconds)

    # Create a context.
    context = mm.Context(system, integrator)
    context.setPositions(pdb.positions)
    
    # Check to make sure energy is finite.
    state = context.getState(getEnergy=True)
    potential = state.getPotentialEnergy() / unit.kilocalories_per_mole
    if numpy.isnan(potential):
        raise Exception("Initial energy for %s is NaN." % pdbcode)

    # Minimize.
    tolerance = 1.0 * unit.kilocalories_per_mole / unit.angstroms
    maxIterations = 50
    mm.LocalEnergyMinimizer.minimize(context, tolerance, maxIterations)

    # Check to make sure energy is finite.
    state = context.getState(getEnergy=True)
    potential = state.getPotentialEnergy() / unit.kilocalories_per_mole
    if numpy.isnan(potential):
        raise Exception("Energy for %s is NaN after minimization." % pdbcode)

    # Simulate.
    nsteps = 500
    integrator.step(nsteps)

    # Check to make sure energy is finite.
    state = context.getState(getEnergy=True)
    potential = state.getPotentialEnergy() / unit.kilocalories_per_mole
    if numpy.isnan(potential):
        raise Exception("Energy for %s is NaN after simulation." % pdbcode)

    del context, integrator

    print "Simulation completed: potential = %.3f kcal/mol" % potential

    return

def test_build_and_simulate():
    # These are tough PDB codes from http://www.umass.edu/microbio/chime/pe_beta/pe/protexpl/badpdbs.htm
    pdbcodes_to_build = ['1AS5', '1CBN', '1DPO', '1IGY', '1HAG', '1IAO', '4CPA', '1QCQ']
    pdbcodes_to_build = ['1VII'] # DEBUG
    pdbcodes_to_simulate = ['1VII'] # should be a subset of pdbcodes_to_build

    # Set up PDB retrieval.
    from Bio.PDB import PDBList
    pdblist = PDBList(server=PDBList.alternative_download_url)

    success = True

    for pdbcode in pdbcodes_to_build:
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
        from simtk import unit
        pH = 7.0
        ionic = 50.0 * unit.millimolar
        box = 10.0 * unit.angstrom
        positiveIon = 'Na+'
        negativeIon = 'Cl-'

        infile = open(input_pdb_filename)
        outfile = open(output_pdb_filename, 'w')

        try:
            from pdbfixer.pdbfixer import PDBFixer, PdbStructure
            from simtk.openmm import app
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
            
            # Test simulating this with OpenMM.
            if pdbcode in pdbcodes_to_simulate:
                simulate(pdbcode, output_pdb_filename)

            # Delete input file.
            os.remove(input_pdb_filename)
            os.remove(output_pdb_filename)

        except Exception as e:
            print str(e)
            success = False

    if not success:
        raise Exception("build test failed on one or more PDB files.")

if __name__ == '__main__':
    test_build_and_simulate()
