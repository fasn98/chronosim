"""Phase 2 gaps found by the scalpel review: functions with no direct test, and F6 checks tighter
than the ones in test_kim_thorne_f6.py. Each tolerance sits above the value measured on 2026-11-07,
quoted next to it.
"""

import mpmath as mp
import pytest

from chronosim import estimates as est
from chronosim import kim_thorne as kt
from chronosim import unruh
from chronosim import wormhole as wh
from chronosim.provenance import Origin

D, B, X = mp.mpf(1000), mp.mpf(10), mp.mpf(250)


def _single_term(xi):
    return -(B / (2 * D)) * (1 + 4 * xi + xi ** 2) / (1 - xi) ** 2 / (6 * mp.pi)


# series_terms: the number of terms of eqs. 66-67 that every numerical stress uses.

@pytest.mark.parametrize("b,xi", [("0.01", "0.5"), ("0.1", "0.9"), ("3", "0.5")])
@pytest.mark.parametrize("tol", [None, "1e-10", "1e-60"])
def test_series_terms_is_the_smallest_N_whose_term_is_below_the_tolerance(b, xi, tol):
    b, xi = mp.mpf(b), mp.mpf(xi)
    N = kt.series_terms(D=1, b=b, xi=xi, tol=tol)
    limit = kt.SERIES_TOL.value if tol is None else mp.mpf(tol)
    ratio = b * xi / 2
    assert ratio ** N <= limit
    assert ratio ** (N - 1) > limit


# unruh_temperature

def test_unruh_temperature_is_one_over_two_pi_xi():
    for xi in ("0.5", "1", "2.5"):
        assert abs(unruh.unruh_temperature(mp.mpf(xi)) * 2 * mp.pi * mp.mpf(xi) - 1) < mp.mpf("1e-75")


# metric_perturbation (eq. 9c) and its relation to the Hawking cutoff

def test_metric_perturbation_is_the_product_of_its_factors():
    q = est.metric_perturbation(mp.mpf("0.1"), 2, mp.mpf("1e-3"), mp.mpf("1e-4"))
    assert abs(q.value / mp.mpf("1e-9") - 1) < mp.mpf("1e-70")
    assert q.origin is Origin.DERIVED


def test_hawking_cutoff_is_where_eq9c_reaches_one():
    dt = est.hawking_cutoff_dt().value
    reached = est.metric_perturbation(1, 0, est.lp_over_D().value, 1 / dt).value
    assert abs(reached - 1) < mp.mpf("1e-12")  # same 1e-12 as the findroot tolerance used in test_estimates.py


# F6, eqs. 54 and 68: the sum over N, which test_kim_thorne_f6.py never exercises.

def test_eq54_with_one_term_is_the_printed_single_term_exactly():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    assert abs(kt.momentum_horizon(mp.mpf(1), D, B, xi, 1).value / _single_term(xi) - 1) < mp.mpf("1e-75")


def test_eq54_is_the_explicit_sum_of_the_single_terms():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    total = sum(
        (B / (2 * D)) ** N * (1 + 4 * xi ** N + xi ** (2 * N)) / (1 - xi ** N) ** 2 for N in (1, 2, 3, 4)
    )
    assert abs(kt.momentum_horizon(mp.mpf(1), D, B, xi, 4).value + total / (6 * mp.pi)) < mp.mpf("1e-75")


def test_eq54_is_inversely_proportional_to_delta_T_and_negative_for_every_order():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    one = kt.momentum_horizon(mp.mpf(1), D, B, xi, 3).value
    assert abs(kt.momentum_horizon(mp.mpf("0.25"), D, B, xi, 3).value / one - 4) < mp.mpf("1e-75")
    assert kt.momentum_horizon(mp.mpf(1), D, B, xi, 3).value < kt.momentum_horizon(mp.mpf(1), D, B, xi, 2).value < 0


def test_eq68_numeric_for_orders_one_to_three_equals_the_eq54_sum():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    numeric = wh.momentum_plane_numeric(mp.mpf(1), X, D, B, "0.5", orders=[1, 2, 3]).value
    assert abs(numeric / kt.momentum_horizon(mp.mpf(1), D, B, xi, 3).value - 1) < mp.mpf("1e-9")  # measured 5.3e-11


def test_eq54_throat_integral_converges_linearly_in_delta_T():
    # The in-repo test stops at delta_T / b = 1e-3 with a bound of 0.01 (measured 0.0091), which a 0.1%
    # error in the coefficient would not cross. Here: measured 0.0787, 0.00907, 0.00102 at 0.1, 0.01, 0.001.
    xi = kt.doppler_xi(mp.mpf("0.5"))
    single = _single_term(xi)
    dev = {}
    for dt in ("0.01", "0.001"):
        v = wh.momentum_throat_numeric(mp.mpf(dt), D, B, "0.5", lambda th: kt.time_Hminus_throat(th, B), orders=[1]).value
        dev[dt] = abs(v / (single / mp.mpf(dt)) - 1)
    assert dev["0.001"] < mp.mpf("1.1e-3")
    assert 5 < dev["0.01"] / dev["0.001"] < 12  # measured 8.87
