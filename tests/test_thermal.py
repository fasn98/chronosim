"""Main arbiter of the factor 1/2: a thermal state in Minkowski space, against Stefan-Boltzmann.

Each tolerance sits a little above the measured error, quoted next to it.
"""

import mpmath as mp
import pytest

from chronosim import thermal, unruh
from chronosim.provenance import Origin

BETA = mp.mpf(1)


@pytest.mark.parametrize("T", ["0.3", "1", "3"])
def test_with_one_half_gives_stefan_boltzmann_and_eq4a_as_printed_twice_it(T):
    T = mp.mpf(T)
    ref = unruh.thermal_density_closed(T)
    out = thermal.thermal_stress(T)
    assert abs(out["with_half"]["rho"].value / ref - 1) < mp.mpf("1e-19")  # measured 3.2e-20
    assert abs(out["eq4a"]["rho"].value / ref - 2) < mp.mpf("1e-18")  # measured 6.4e-20


@pytest.mark.parametrize("name", ["eq4a", "with_half"])
def test_pressure_is_a_third_of_the_density_and_the_trace_vanishes(name):
    d = thermal.thermal_stress(1)[name]
    rho = d["rho"].value
    for p in d["p"]:
        assert abs(3 * p.value / rho - 1) < mp.mpf("1e-19")  # measured 2.5e-20
    assert abs(d["trace"].value / rho) < mp.mpf("1e-21")  # measured 8.4e-23
    assert d["max_offdiag"].value / abs(rho) < mp.mpf("1e-30")  # measured 0
    assert d["rho"].origin is Origin.COMPUTED


def test_density_has_dimension_length_to_minus_4():
    # The split scales with beta, so the ratio is exact up to rounding. That floor is about
    # 1e-80 / (h / beta)^2 = 1e-58 for the second differences. Measured 2.5e-61.
    assert abs(thermal.thermal_stress(3)["with_half"]["rho"].value / thermal.thermal_stress(1)["with_half"]["rho"].value / 81 - 1) < mp.mpf("1e-59")


def test_terms_decay_as_n_to_minus_4_and_the_tail_is_needed():
    d = thermal.thermal_stress(1)["with_half"]
    assert d["decay_spread"].value < mp.mpf("1e-21")  # measured 1.8e-22
    assert mp.mpf("1e-5") < d["tail_fraction"].value < mp.mpf("3e-5")  # measured 2.1e-5


def test_result_does_not_depend_on_the_point_or_on_the_number_of_explicit_terms():
    a = thermal.thermal_stress(1)["with_half"]["rho"].value
    moved = thermal.thermal_stress(1, X=["0.3", "0.2", "-0.5", "1.1"])["with_half"]["rho"].value
    assert abs(moved / a - 1) < mp.mpf("1e-60")  # measured 0
    assert abs(thermal.thermal_stress(1, n_terms=12)["with_half"]["rho"].value / a - 1) < mp.mpf("1e-25")  # measured 1.6e-26


@pytest.mark.parametrize("t,r,tol", [("0.1", "0.2", "1e-9"), ("0.3", "0.5", "1e-8"), ("0.7", "1.5", "1e-7")])
def test_image_sum_matches_the_mode_integral(t, r, tol):
    # Measured 2.3e-10, 2.1e-9 and 2.5e-8 with 400 pairs. The difference is the cut of the image sum.
    t, r, N = mp.mpf(t), mp.mpf(r), 400
    partial = thermal.hadamard_regularized([t, r, 0, 0], [0, 0, 0, 0], BETA, N)
    tail = (mp.zeta(2) - sum(mp.mpf(n) ** -2 for n in range(1, N + 1))) / (mp.pi ** 2 * BETA ** 2)
    assert abs((partial + tail) / thermal.hadamard_mode_integral(t, r, BETA).value - 1) < mp.mpf(tol)


def test_coincidence_value_of_the_image_sum():
    # (1 / pi^2) int k n_k dk = T^2 / 6, and the pair sum at coincidence is zeta(2) / pi^2 for beta = 1.
    integral = mp.quad(lambda k: k / mp.expm1(k), [0, 1, 10, mp.inf]) / mp.pi ** 2
    assert abs(integral - mp.mpf(1) / 6) < mp.mpf("1e-14")
    assert abs(mp.zeta(2) / mp.pi ** 2 - integral) < mp.mpf("1e-14")
