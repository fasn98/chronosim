"""Inputs outside the domain raise ValueError instead of returning a silent wrong number (reviewer I4)."""

import mpmath as mp
import pytest

from chronosim import kim_thorne as kt
from chronosim import thermal, unruh
from chronosim.provenance import Origin

D, B, X = mp.mpf(1), mp.mpf("0.01"), mp.mpf("0.25")
XI = kt.doppler_xi(mp.mpf("0.5"))
P = [mp.mpf("0.2"), X, mp.mpf("0.01"), mp.mpf(0)]
Q = [mp.mpf("0.19"), mp.mpf("0.26"), mp.mpf("0.02"), mp.mpf(0)]


@pytest.mark.parametrize("beta", ["0", "1", "-0.5", "1.5"])
def test_doppler_xi_requires_0_lt_beta_lt_1(beta):
    with pytest.raises(ValueError):
        kt.doppler_xi(mp.mpf(beta))


@pytest.mark.parametrize("b,xi", [("2", "1"), ("4", "0.5"), ("0", "0.5")])
def test_series_terms_requires_a_ratio_in_0_1(b, xi):
    # b xi / 2D = 1, 1 and 0: a non-convergent or degenerate series.
    with pytest.raises(ValueError):
        kt.series_terms(D, mp.mpf(b), mp.mpf(xi))


@pytest.mark.parametrize("call", [
    lambda: kt.hadamard_between(P, Q, D, B, XI, 0),
    lambda: kt.hadamard_throat(P, Q, D, B, XI, 0),
    lambda: kt.stress_between_kk(mp.mpf("0.2"), X, mp.mpf(0), D, B, XI, 0),
    lambda: kt.stress_between_kk(mp.mpf("0.2"), X, mp.mpf(0), D, B, XI, 3, []),
    lambda: kt.stress_between_kk(mp.mpf("0.2"), X, mp.mpf(0), D, B, XI, 3, [0]),
    lambda: kt.stress_throat(mp.mpf("-0.001"), mp.mpf("0.7"), mp.mpf(0), D, B, XI, 0, "transverse"),
    lambda: kt.stress_throat(mp.mpf("-0.001"), mp.mpf("0.7"), mp.mpf(0), D, B, XI, 3, "transverse", []),
    lambda: kt.klinkhammer_kk_on_axis(mp.mpf("1e-3"), X, D, B, XI, 0),
    lambda: kt.momentum_horizon(mp.mpf(1), D, B, XI, 0),
    lambda: kt.zeta(0, D, XI),
    lambda: kt.zeta_between_sum(0, X, D, XI),
    lambda: thermal.hadamard_regularized([0, 1, 0, 0], [0, 0, 0, 0], mp.mpf(1), 0),
    lambda: thermal.thermal_stress(1, n_terms=1),
], ids=["hadamard_between", "hadamard_throat", "stress_between_kk", "kk_empty_orders", "kk_order_0",
        "stress_throat", "throat_empty_orders", "klinkhammer", "momentum_horizon", "zeta", "zeta_between_sum",
        "thermal_images", "thermal_stress"])
def test_empty_sums_raise(call):
    with pytest.raises(ValueError):
        call()


def test_a_sum_over_1_to_n_terms_is_computed_even_with_one_term():
    # Reviewer note: series_terms may return 1, and that single term is still a truncated series.
    q = kt.stress_between_kk(mp.mpf("0.2"), X, mp.mpf(0), D, B, XI, 1)
    assert q.origin is Origin.COMPUTED
    assert kt.stress_between_kk(mp.mpf("0.2"), X, mp.mpf(0), D, B, XI, 1, [1]).origin is Origin.DERIVED


def test_unruh_reference_is_heuristic_with_its_caveat():
    q = unruh.rindler_density(mp.mpf(1))
    assert q.origin is Origin.HEURISTIC
    assert "article not read" in q.source


@pytest.mark.parametrize("n_terms", [1, 3])
def test_momentum_horizon_and_stress_between_kk_follow_the_same_origin_rule(n_terms):
    # A sum over 1..n_terms is COMPUTED even with one term; one explicit term through orders is DERIVED.
    sum_kk = kt.stress_between_kk(mp.mpf("0.2"), X, mp.mpf(0), D, B, XI, n_terms)
    sum_p = kt.momentum_horizon(mp.mpf(1), D, B, XI, n_terms)
    assert sum_kk.origin is sum_p.origin is Origin.COMPUTED
    one_kk = kt.stress_between_kk(mp.mpf("0.2"), X, mp.mpf(0), D, B, XI, n_terms, [1])
    one_p = kt.momentum_horizon(mp.mpf(1), D, B, XI, n_terms, [1])
    assert one_kk.origin is one_p.origin is Origin.DERIVED
