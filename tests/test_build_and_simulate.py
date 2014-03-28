
import pdbfixer
import simtk.openmm 
import Bio.PDB

import os
import os.path
import sys
import numpy

from threading import Timer

#from a solution on stackoverflow
class Watchdog:
    def __init__(self, timeout, userHandler=None):  # timeout in seconds
        self.timeout = timeout
        self.handler = userHandler if userHandler is not None else self.defaultHandler
        self.timer = Timer(self.timeout, self.handler)

    def reset(self):
        self.timer.cancel()
        self.timer = Timer(self.timeout, self.handler)

    def stop(self):
        self.timer.cancel()

    def defaultHandler(self):
        raise self

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
    pdbcodes_to_simulate = ['1AS5', '1CBN', '1DPO', '1IGY', '1HAG', '1IAO', '4CPA', '1QCQ']

    # Set up PDB retrieval.
    from Bio.PDB import PDBList
    pdblist = PDBList(server=PDBList.alternative_download_url)

    success = True

    for pdbcode in pdbcodes_to_build:
        print pdbcode

        try:
            # Remove file if it exists already.
            input_pdb_filename = 'pdb' + pdbcode + '.ent'
            if os.path.exists(input_pdb_filename):  
                os.remove(input_pdb_filename)                

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
        
        success = True

        timeout_seconds = 0.1
        watchdog = Watchdog(timeout_seconds)
        try:        
            from pdbfixer.pdbfixer import PDBFixer, PdbStructure
            from simtk.openmm import app
            stage = "Creating PDBFixer..."
            fixer = PDBFixer(PdbStructure(infile))
            stage = "Finding missing residues..."
            fixer.findMissingResidues()
            stage = "Finding nonstandard residues..."
            fixer.findNonstandardResidues()
            stage = "Replacing nonstandard residues..."
            fixer.replaceNonstandardResidues()
            stage = "Finding missing atoms..."
            fixer.findMissingAtoms()
            stage = "Adding missing atoms..."
            fixer.addMissingAtoms()
            stage = "Removing heterogens..."
            fixer.removeHeterogens(False)
            stage = "Adding missing hydrogens..."
            fixer.addMissingHydrogens(pH)
            stage = "Writing PDB file..."
            app.PDBFile.writeFile(fixer.topology, fixer.positions, outfile)
            stage = "Done."
            infile.close()
            outfile.close()

        except Watchdog:
            print "timed out fixing PDB %s" % pdbcode
            success = False

        except Exception as e:
            print str(e)
            success = False
        
        watchdog.stop()
        del watchdog
                    
        # Test simulating this with OpenMM.
        if pdbcode in pdbcodes_to_simulate:
            watchdog = Watchdog(timeout_seconds)
            try:
                simulate(pdbcode, output_pdb_filename)
                
            except Watchdog:
                print "PDB code %s timed out in stage '%s'." % (pdbcode, stage)
                success = False

            except Exception as e:
                print str(e)
                success = False
        
            watchdog.stop()
            del watchdog

        # Clean up.
        os.remove(input_pdb_filename)
        os.remove(output_pdb_filename)

    if not success:
        raise Exception("build test failed on one or more PDB files.")

if __name__ == '__main__':
    test_build_and_simulate()
