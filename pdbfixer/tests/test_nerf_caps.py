"""Tests for NeRF-based ACE/NME cap placement."""

from io import StringIO

import numpy as np
import pytest
from openmm.app import ForceField
from pdbfixer import PDBFixer
from pdbfixer.pdbfixer import _placeAtomNeRF, _buildAceCap, _buildNmeCap


def _dist(a, b):
    return np.linalg.norm(a - b)


def _angle(a, b, c):
    ba = a - b
    bc = c - b
    cos_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    return np.arccos(np.clip(cos_angle, -1, 1))


def _dihedral(a, b, c, d):
    b1 = b - a
    b2 = c - b
    b3 = d - c
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m = np.cross(n1, b2 / np.linalg.norm(b2))
    return np.arctan2(np.dot(m, n2), np.dot(n1, n2))


# Reference backbone atoms for a typical residue (in nm).
REF_N  = np.array([0.0, 0.0, 0.0])
REF_CA = np.array([0.1458, 0.0, 0.0])
REF_C  = np.array([0.2058, 0.1230, 0.0])


def _rotation_matrix(axis, angle):
    """Rodrigues' rotation formula."""
    axis = axis / np.linalg.norm(axis)
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0],
    ])
    return np.identity(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K


@pytest.fixture(params=[
    (np.identity(3), np.zeros(3)),
    (_rotation_matrix(np.array([1.0, 2.0, 3.0]), 1.23), np.array([0.5, -0.3, 0.7])),
], ids=["z=0 plane", "general 3D"])
def backbone(request):
    """Yield (N, CA, C) reference backbone, optionally rotated+translated."""
    rot, shift = request.param
    return rot @ REF_N + shift, rot @ REF_CA + shift, rot @ REF_C + shift


class TestPlaceAtomNeRF:
    # Use non-collinear reference points (collinear points make cross product zero).
    A1 = np.array([0.0, 0.1, 0.0])
    A2 = np.array([0.0, 0.0, 0.0])
    A3 = np.array([0.15, 0.0, 0.0])

    def test_bond_length_and_angle(self):
        bond_length = 0.1335
        bond_angle = 2.035
        x = _placeAtomNeRF(self.A1, self.A2, self.A3, bond_length, bond_angle, np.pi / 3)
        assert abs(_dist(self.A3, x) - bond_length) < 1e-6
        assert abs(_angle(self.A2, self.A3, x) - bond_angle) < 1e-6

    def test_collinear_raises(self):
        a1 = np.array([0.0, 0.0, 0.0])
        a2 = np.array([0.15, 0.0, 0.0])
        a3 = np.array([0.30, 0.0, 0.0])
        with pytest.raises(ValueError, match="collinear"):
            _placeAtomNeRF(a1, a2, a3, 0.15, 2.0, 1.0)

    def test_dihedral_pi(self):
        """Test with pi dihedral, which is sign-invariant and matches cap usage."""
        x = _placeAtomNeRF(self.A1, self.A2, self.A3, 0.15, 2.0, np.pi)
        measured = _dihedral(self.A1, self.A2, self.A3, x)
        assert abs(abs(measured) - np.pi) < 1e-5


class TestBuildAceCap:
    def test_ace_bond_lengths(self, backbone):
        n, ca, c = backbone
        caps = _buildAceCap(n, ca, c)
        assert abs(_dist(caps['C'], n) - 0.1329) < 1e-4
        assert abs(_dist(caps['O'], caps['C']) - 0.1231) < 1e-4
        assert abs(_dist(caps['CH3'], caps['C']) - 0.1525) < 1e-4

    def test_ace_omega_trans(self, backbone):
        n, ca, c = backbone
        caps = _buildAceCap(n, ca, c)
        omega = _dihedral(caps['CH3'], caps['C'], n, ca)
        assert abs(abs(omega) - np.pi) < 0.05


class TestBuildNmeCap:
    def test_nme_bond_lengths(self, backbone):
        n, ca, c = backbone
        caps = _buildNmeCap(n, ca, c)
        assert abs(_dist(caps['N'], c) - 0.1329) < 1e-4
        assert abs(_dist(caps['C'], caps['N']) - 0.1458) < 1e-4

    def test_nme_omega_trans(self, backbone):
        n, ca, c = backbone
        caps = _buildNmeCap(n, ca, c)
        omega = _dihedral(ca, c, caps['N'], caps['C'])
        assert abs(abs(omega) - np.pi) < 0.05


class TestCapIntegration:
    def test_ace_nme_forcefield_compatible(self):
        """Full PDBFixer pipeline: add caps and verify ForceField accepts them."""
        pdb_text = """\
SEQRES   1 A    4  ACE ALA ALA NME
ATOM      1  N   ALA A   2       1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   2       2.009   1.420   0.000  1.00  0.00           C
ATOM      3  C   ALA A   2       3.530   1.411   0.000  1.00  0.00           C
ATOM      4  O   ALA A   2       4.135   0.348   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   2       1.499   2.156   1.231  1.00  0.00           C
ATOM      6  N   ALA A   3       4.090   2.600   0.000  1.00  0.00           N
ATOM      7  CA  ALA A   3       5.540   2.700   0.000  1.00  0.00           C
ATOM      8  C   ALA A   3       6.100   4.100   0.000  1.00  0.00           C
ATOM      9  O   ALA A   3       5.400   5.100   0.000  1.00  0.00           O
ATOM     10  CB  ALA A   3       6.050   1.964   1.231  1.00  0.00           C
END
"""
        fixer = PDBFixer(pdbfile=StringIO(pdb_text))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)

        assert fixer.topology.getNumAtoms() == 32

        # If geometry is bad, ForceField will raise an exception.
        ff = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        ff.createSystem(fixer.topology)
