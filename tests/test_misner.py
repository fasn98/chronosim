import mpmath as mp
import pytest

from chronosim.misner import (
    closed_form_stress,
    hadamard_image,
    hadamard_minkowski,
    self_consistent_b,
    stress_prefactor,
    to_rindler,
)
from chronosim.provenance import Origin


def test_rindler_coordinates_round_trip():
    xi, eta = mp.mpf("1.7"), mp.mpf("0.3")
    t, x = xi * mp.sinh(eta), xi * mp.cosh(eta)
    xi2, eta2 = to_rindler(t, x)
    assert abs(xi2 - xi) < mp.mpf("1e-70") and abs(eta2 - eta) < mp.mpf("1e-70")


def test_points_outside_region_r_are_rejected():
    for t, x in [(1, 1), (2, 1), (0, -1)]:
        with pytest.raises(ValueError):
            to_rindler(t, x)


def test_rindler_term_has_the_minkowski_short_distance_singularity():
    # Hadamard form: G_0 / G_M -> 1 as the points merge
    X = [0, 1, 0, 0]
    for d in ["1e-4", "1e-8", "1e-16"]:
        Xp = [0, 1, 0, mp.mpf(d)]
        r = hadamard_image(0, X, Xp, 1) / hadamard_minkowski(X, Xp)
        assert abs(r - 1) < 10 * mp.mpf(d)


def test_images_are_symmetric_under_n_to_minus_n_at_equal_time():
    X, Xp = [0, 1, 0, 0], [0, 1, 0, mp.mpf("0.1")]
    assert abs(hadamard_image(3, X, Xp, 2) - hadamard_image(-3, X, Xp, 2)) < mp.mpf("1e-70")


def test_closed_form_is_traceless_and_derived():
    T = closed_form_stress(1, 1)
    assert -T["rho"].value + T["p_xi"].value + T["p_y"].value + T["p_z"].value == 0
    assert all(q.origin is Origin.DERIVED for q in T.values())
    assert all("eq. 12" in q.source for q in T.values())


def test_closed_form_vanishes_at_b_equals_2pi():
    b = self_consistent_b()
    assert b.origin is Origin.CITED
    assert abs(stress_prefactor(1, b.value)) < mp.mpf("1e-75")


def test_closed_form_rindler_limit_and_horizon_scaling():
    # b -> infinity: Rindler vacuum, rho = -3 / (1440 pi^2 xi^4)
    assert abs(closed_form_stress(1, mp.mpf("1e30"))["rho"].value + 3 / (1440 * mp.pi ** 2)) < mp.mpf("1e-60")
    # xi^-4: halving the distance to the horizon multiplies the energy by 16
    assert abs(stress_prefactor("0.5", 3) / stress_prefactor(1, 3) - 16) < mp.mpf("1e-70")
