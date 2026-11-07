"""Arbiter for the factor of 2: Kim & Thorne eqs. 3, 4 and 20 on the Rindler vacuum.

Each tolerance sits a little above the error measured at delta = 1e-10, h = 1e-11,
dps = 80, and the measured value is quoted next to it.
"""

import mpmath as mp
import pytest

from chronosim.kim_thorne import hadamard_flat
from chronosim.misner import hadamard_image, hadamard_minkowski
from chronosim.pointsplit import rindler_part
from chronosim.wormhole import eq4_stress, half_square_separation, regularized_with_eq3

EXPECTED_PHASE1 = (-3, -1, -1, -1)  # units of 1 / (1440 pi^2 xi^4), validated in phase 1


def _rindler_eq4(xi):
    F = regularized_with_eq3(lambda P, Q: hadamard_image(0, P, Q, 1))
    return eq4_stress(F, [0, xi, 0, 0])


def test_operator_coefficients_on_a_known_function():
    # F = t t': 2/3 - 1/6 * g_tt * (-1) = 1/2 for T_tt, -1/6 * g_xx * (-1) = 1/6 for T_xx, no factor 1/2.
    T = eq4_stress(lambda P, Q: P[0] * Q[0], [0, 1, 0, 0])
    assert abs(T[0][0] - mp.mpf(1) / 2) < mp.mpf("1e-30")
    assert abs(T[1][1] - mp.mpf(1) / 6) < mp.mpf("1e-30")


def test_the_two_articles_normalize_the_flat_hadamard_function_identically():
    X, Q = [0, mp.mpf("1.3"), 0, 0], [0, mp.mpf("1.3"), 0, mp.mpf("1e-6")]
    sigma0 = half_square_separation(X, Q)
    assert abs(hadamard_flat(sigma0) / hadamard_minkowski(X, Q) - 1) < mp.mpf("1e-75")


def test_rindler_function_has_the_kim_thorne_singularity():
    # G_R 4 pi^2 sigma_0 -> 1, the N = 0 term of eq. 20. Measured 1 - 4.9e-14 at a split of 1e-6.
    X, Q = [0, mp.mpf("1.3"), 0, 0], [0, mp.mpf("1.3"), 0, mp.mpf("1e-6")]
    assert abs(hadamard_image(0, X, Q, 1) / hadamard_flat(half_square_separation(X, Q)) - 1) < mp.mpf("1e-13")


@pytest.mark.parametrize("xi", ["0.5", "1.3", "2"])
def test_eq4a_gives_twice_the_validated_rindler_stress(xi):
    xi = mp.mpf(xi)
    scale = 1 / (1440 * mp.pi ** 2 * xi ** 4)
    T = _rindler_eq4(xi)
    R = rindler_part([0, xi, 0, 0])
    for i in range(4):
        assert abs(R[i][i] / scale - EXPECTED_PHASE1[i]) < mp.mpf("1e-15")  # the phase 1 tolerance
        assert abs(T[i][i] / R[i][i] - 2) < mp.mpf("1e-18")  # measured at most 7.3e-19
