"""Arbiter independent of both articles: the Rindler vacuum against the Unruh value."""

import mpmath as mp
import pytest

from chronosim import unruh
from chronosim.misner import hadamard_image
from chronosim.pointsplit import rindler_part
from chronosim.provenance import Origin
from chronosim.wormhole import eq4_stress, regularized_with_eq3


def test_thermal_density_against_the_bose_integral():
    T = mp.mpf("0.3")
    q = unruh.thermal_density_integral(T)
    assert q.origin is Origin.COMPUTED
    assert abs(q.value / unruh.thermal_density_closed(T) - 1) < mp.mpf("1e-79")  # measured 2.1e-81


def test_unruh_value_is_the_stated_rational_multiple():
    xi = mp.mpf("1.3")
    assert abs(unruh.rindler_density(xi).value * 1440 * mp.pi ** 2 * xi ** 4 + 3) < mp.mpf("1e-75")


def test_unruh_value_has_dimension_length_to_minus_4():
    assert abs(unruh.rindler_density(3).value / unruh.rindler_density(1).value * 81 - 1) < mp.mpf("1e-75")


@pytest.mark.parametrize("xi", ["0.5", "1.3", "2"])
def test_operator_with_one_half_agrees_and_eq4a_as_printed_is_twice(xi):
    xi = mp.mpf(xi)
    ref = unruh.rindler_density(xi).value
    F = regularized_with_eq3(lambda P, Q: hadamard_image(0, P, Q, 1))
    assert abs(rindler_part([0, xi, 0, 0])[0][0] / ref - 1) < mp.mpf("1e-18")  # measured at most 5.8e-20
    assert abs(eq4_stress(F, [0, xi, 0, 0])[0][0] / ref - 2) < mp.mpf("1e-18")  # measured at most 1.2e-19
