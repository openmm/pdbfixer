
from __future__ import print_function
from pdbfixer.pdbfixer import PDBFixer
from openmm import app

import os
import os.path
import sys
import numpy
import tempfile

from threading import Timer

#from a solution on stackoverflow
class Watchdog(BaseException):
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
    from openmm import app
    import openmm.openmm as mm
    from openmm import unit
    from sys import stdout

    # Load the PDB file.
    pdb = app.PDBFile(pdb_filename)
    
    # Set up implicit solvent forcefield.
    #forcefield = app.ForceField('amber99sbildn.xml')
    forcefield = app.ForceField('amber10.xml')
    
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

    print("Simulation completed: potential = %.3f kcal/mol" % potential)

    return

def test_build_and_simulate():
    # DEBUG: These are tough PDB codes from http://www.umass.edu/microbio/chime/pe_beta/pe/protexpl/badpdbs.htm
    pdbcodes_to_build = ['1AS5', '1CBN', '1DPO', '1IGY', '1HAG', '1IAO', '4CPA', '1QCQ']

    # DEBUG: Small test cases.
    pdbcodes_to_build = ['110D', '116D', '117D', '118D', '134D', '135D', '136D', '138D', '143D', '148D', '151D', '152D', '159D', '177D', '17RA', '183D', '184D', '186D', '187D', '188D', '189D', '1A11', '1A13', '1A1P', '1A3P', '1A51', '1A60', '1A83', '1A9L', '1AAF', '1AB1', '1ABZ', '1AC7', '1ACW', '1AD7', '1ADX', '1AFP', '1AFT', '1AFX', '1AG7', '1AGG', '1AGL', '1AGT', '1AHL', '1AIE', '1AJ1', '1AJF', '1AJJ', '1AJU', '1AKG', '1AKX', '1AL1', '1ALE', '1ALF', '1ALG', '1AM0', '1AMB', '1AMC', '1AML', '1ANP', '1ANR', '1ANS', '1AO9', '1AOO', '1APF', '1APO', '1APQ', '1AQG', '1AQO', '1AQQ', '1AQR', '1AQS', '1ARD', '1ARE', '1ARF', '1ARJ', '1ARK', '1AS5', '1AT4', '1ATO', '1ATV', '1ATW', '1ATX', '1AV3', '1AW4', '1AW6', '1AWY', '1AXH', '1AY3', '1AYJ', '1AZ6', '1AZH', '1AZJ', '1AZK', '1B03', '1B0Q', '1B13', '1B1V', '1B2J', '1B36', '1B45', '1B4G', '1B4I', '1B4Y', '1B5N', '1B8W', '1B9G', '1B9P', '1B9Q', '1B9U', '1BA4', '1BA5', '1BA6', '1BAH', '1BAL', '1BBA', '1BBG', '1BBL', '1BBO', '1BCV', '1BD1', '1BDC', '1BDD', '1BDE', '1BDK', '1BDS', '1BE7', '1BEI', '1BF0', '1BF9', '1BFW', '1BFY', '1BFZ', '1BGK', '1BGZ', '1BH0', '1BH1', '1BH4', '1BH7', '1BHI', '1BHP', '1BIG', '1BJB', '1BJC', '1BJH', '1BK2', '1BK8', '1BKT', '1BKU', '1BL1', '1BM4', '1BMX', '1BN0', '1BNB', '1BNX', '1BOE', '1BOR', '1BPI', '1BPT', '1BQ8', '1BQ9', '1BQF', '1BRF', '1BRV', '1BRZ', '1BTI', '1BTQ', '1BTR', '1BTS', '1BTT', '1BUB', '1BUS', '1BVJ', '1BW6', '1BWX', '1BX7', '1BX8', '1BY0', '1BY6', '1BYJ', '1BYV', '1BYY', '1BZ2', '1BZ3', '1BZB', '1BZG', '1BZK', '1BZT', '1BZU', '1C0O', '1C26', '1C2U', '1C32', '1C34', '1C35', '1C38', '1C49', '1C4B', '1C4E', '1C4S', '1C55', '1C56', '1C6W', '1C98', '1C9A', '1C9Z', '1CAA', '1CAD', '1CAP', '1CB3', '1CB9', '1CBH', '1CBN', '1CCF', '1CCM', '1CCN', '1CCQ', '1CCV', '1CE3', '1CE4', '1CEK', '1CEU', '1CFG', '1CFH', '1CFI', '1CHL', '1CHV', '1CIX', '1CKW', '1CKX', '1CKY', '1CKZ', '1CL4', '1CLF', '1CMR', '1CNL', '1CNN', '1CNR', '1CO4', '1COI', '1CQ0', '1CQ5', '1CQL', '1CQU', '1CR8', '1CRE', '1CRF', '1CRN', '1CS9', '1CSA', '1CT6', '1CTI', '1CV9', '1CVQ', '1CW5', '1CW6', '1CW8', '1CWX', '1CWZ', '1CXN', '1CXO', '1CXR', '1CXW', '1CYA', '1CYB', '1CZ6', '1D0R', '1D0T', '1D0U', '1D0W', '1D10', '1D11', '1D12', '1D13', '1D14', '1D15', '1D16', '1D17', '1D1E', '1D1F', '1D1H', '1D26', '1D2D', '1D2J', '1D2L', '1D33', '1D35', '1D36', '1D37', '1D38', '1D54', '1D58', '1D5Q', '1D61', '1D62', '1D67', '1D6B', '1D6X', '1D78', '1D79', '1D7N', '1D7T', '1D7Z', '1D82', '1D8G', '1D93', '1D9J', '1D9L', '1D9M', '1D9O', '1D9P', '1DA0', '1DA9', '1DB6', '1DEC', '1DEM', '1DEN', '1DEP', '1DF6', '1DFE', '1DFS', '1DFT', '1DFW', '1DFY', '1DFZ']

    # impossible cases
    pdbcodes_to_build = [
        '1AO9', # contains residue DOP, which is not resolved in the ATOM records and does not appear to have a machine-readable definition anywhere
        ]

    # DEBUG: A few small test cases.
    pdbcodes_to_build = ['110D', '116D', '117D', '118D', '134D', '135D', '136D', '138D', '143D', '148D', '151D', '152D', '177D', '17RA', '183D', '184D', '186D', '187D', '188D', '189D', '1A11', '1A13', '1A1P', '1A3P', '1A51', '1A60', '1A83', '1A9L', '1AAF', '1AB1', '1ABZ', '1AC7', '1ACW', '1AD7', '1ADX', '1AFP', '1AFT', '1AFX', '1AG7', '1AGG', '1AGL', '1AGT', '1AHL', '1AIE', '1AJ1', '1AJF', '1AJJ', '1AKG', '1AL1', '1ALE', '1ALF', '1ALG', '1AM0', '1AMB', '1AMC', '1AML', '1ANP', '1ANR', '1ANS', '1AOO', '1BH7', '1BX8', '1CEK']

    # Don't simulate any.
    pdbcodes_to_simulate = []

    # Keep track of list of failures.
    failures = list()
        
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    for pdbcode in pdbcodes_to_build:
        print("------------------------------------------------")
        print(pdbcode)

        pH = 7.0

        outfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
        output_pdb_filename = outfile.name

        timeout_seconds = 30
        watchdog = Watchdog(timeout_seconds)
        build_successful = False
        try:
            stage = "Creating PDBFixer..."
            fixer = PDBFixer(pdbid=pdbcode)
            stage = "Deleting hydrogens..."
            if pdbcode in ['135D', '136D', '177D', '1A83', '1AGG', '1AJ1']:
                # These input files include extra hydrogens that aren't supported by the force field.
                # To avoid problems, delete all pre-existing hydrogens.
                modeller = app.Modeller(fixer.topology, fixer.positions)
                modeller.delete([a for a in fixer.topology.atoms() if a.element == app.element.hydrogen])
                fixer.topology = modeller.topology
                fixer.positions = modeller.positions
            stage = "Finding missing residues..."
            fixer.findMissingResidues()
            stage = "Finding nonstandard residues..."
            fixer.findNonstandardResidues()
            stage = "Replacing nonstandard residues..."
            fixer.replaceNonstandardResidues()
            stage = "Removing heterogens..."
            fixer.removeHeterogens(False)
            stage = "Finding missing atoms..."
            fixer.findMissingAtoms()
            stage = "Adding missing atoms..."
            fixer.addMissingAtoms()
            stage = "Adding missing hydrogens..."
            fixer.addMissingHydrogens(pH)
            stage = "Writing PDB file..."
            app.PDBFile.writeFile(fixer.topology, fixer.positions, outfile)
            stage = "Create System..."
            forcefield.createSystem(fixer.topology)
            stage = "Done."
            outfile.close()
            build_successful = True

        except Watchdog:
            message = "timed out in stage %s" % stage
            print(message)
            failures.append((pdbcode, Exception(message)))

        except Exception as e:
            print("EXCEPTION DURING BUILD")
            #import traceback
            #print traceback.print_exc()
            print(str(e))
            failures.append((pdbcode, e))
        
        watchdog.stop()
        del watchdog
                    
        # Test simulating this with OpenMM.
        if (pdbcode in pdbcodes_to_simulate) and (build_successful):
            watchdog = Watchdog(timeout_seconds)
            try:
                simulate(pdbcode, output_pdb_filename)
                
            except Watchdog:
                message = "timed out in simulation"
                print(message)
                failures.append((pdbcode, Exception(message)))

            except Exception as e:
                print("EXCEPTION DURING SIMULATE")
                #import traceback
                #print traceback.print_exc()
                print(str(e))
                failures.append((pdbcode, e))
        
            watchdog.stop()
            del watchdog

        # Clean up.
        os.remove(output_pdb_filename)

    print("------------------------------------------------")

    if len(failures) != 0:
        print("")
        print("SUMMARY OF FAILURES:")
        print("")
        for failure in failures:
            (pdbcode, exception) = failure
            print("%6s : %s" % (pdbcode, str(exception)))
        print("")

        raise Exception("Build test failed on one or more PDB files.")
    
    else:
        print("All tests succeeded.")

if __name__ == '__main__':
    test_build_and_simulate()
