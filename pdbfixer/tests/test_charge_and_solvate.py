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
    ],
)
def test_charge_and_solvate(pdbCode, soluteCharge):
    """
    Test that downloadCharges and addSolvent successfully neutralise the system

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

    chains_to_cap = {chain for chain, resi in fixer.missingResidues}
    for chainidx in chains_to_cap:
        chain = [*fixer.topology.chains()][chainidx]
        last_resi = len([*chain.residues()])
        # Capping with GLY because ACE/NME currently breaks addMissingHydrogens
        # Adding a cap keeps the protein compact and the addSolvent call quick
        fixer.missingResidues[chainidx, 0] = ["GLY"]
        fixer.missingResidues[chainidx, last_resi] = ["GLY"]
        # fixer.missingResidues[chainidx, 0] = ['ACE']
        # fixer.missingResidues[chainidx, last_resi] = ['NME']

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    fixer.addMissingHydrogens(pH=7.4)
    fixer.downloadCharges()
    fixer.addSolvent(
        padding=2.0 * openmm.unit.nanometer,
        ionicStrength=0.1 * openmm.unit.molar,
        boxShape="dodecahedron",
    )

    numCl = sum(1 for res in fixer.topology.residues() if res.name.lower() == "cl")
    numNa = sum(1 for res in fixer.topology.residues() if res.name.lower() == "na")
    assert soluteCharge == numCl - numNa
