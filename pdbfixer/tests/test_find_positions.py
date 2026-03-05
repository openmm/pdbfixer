from openmm import Vec3
from pytest import approx

from pdbfixer.pdbfixer import _findOXTPosition

ATOM_POSITIONS_ALA = {
    "N": Vec3(-0.966, 0.493, 1.500),
    "CA": Vec3(0.257, 0.418, 0.692),
    "C": Vec3(-0.094, 0.017, -0.716),
    "O": Vec3(-1.056, -0.682, -0.923),
    "CB": Vec3(1.204, -0.620, 1.296),
    "OXT": Vec3(0.586, 0.394, -1.639),
}


def test_findOXTPosition():
    pos_oxt = _findOXTPosition(ATOM_POSITIONS_ALA)
    assert pos_oxt == approx(ATOM_POSITIONS_ALA["OXT"], abs=1e-3)
