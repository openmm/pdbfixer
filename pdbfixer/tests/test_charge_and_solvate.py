import pytest
from pdbfixer import PDBFixer
import openmm.unit


@pytest.mark.parametrize(
    "pdbCode,soluteCharge",
    [
        ("1PO0", -21),
        ("1A11", 1),
        ("110D", -5),
        ("1CNR", 0),
        ("1ESD", -21),
        ("25c8", -2),
    ],
)
def test_charge_and_solvate(pdbCode, soluteCharge):
    """
    Test that addSolvent successfully neutralises systems

    Parameters
    ----------
    pdbCode : str
        The PDB ID to test
    soluteCharge : int
        The formal charge of the solute - should equal the number of chloride
        ions minus the number of sodium ions in the neutralised system. Note
        that this may include ions from the original PDB entry.
    """
    fixer = PDBFixer(pdbid=pdbCode)
    fixer.findMissingResidues()

    chainLengths = [len([*chain.residues()]) for chain in fixer.topology.chains()]
    for chainidx, residx in list(fixer.missingResidues):
        if residx == 0:
            fixer.missingResidues[chainidx, residx] = ["GLY"]
        elif residx == chainLengths[chainidx]:
            fixer.missingResidues[chainidx, residx] = ["GLY"]

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    fixer.addSolvent(
        padding=2.0 * openmm.unit.nanometer,
        ionicStrength=0.1 * openmm.unit.molar,
        boxShape="dodecahedron",
    )

    numCl = sum(1 for res in fixer.topology.residues() if res.name.lower() == "cl")
    numNa = sum(1 for res in fixer.topology.residues() if res.name.lower() == "na")
    assert soluteCharge == numCl - numNa
