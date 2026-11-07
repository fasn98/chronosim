"""F3: Kim & Thorne eqs. 61 and 66 -> eq. 67, without using eq. 67 in the numerics.

Each tolerance sits a little above the measured value, quoted next to it.
"""

import mpmath as mp
import pytest

from chronosim import kim_thorne as kt
from chronosim import wormhole as wh
from chronosim.provenance import Origin

D = mp.mpf(1)
COMPONENTS = ("TT", "TX", "XX")


def _axis(beta, b, X, delta, **kw):
    b, X, delta = mp.mpf(b), mp.mpf(X), mp.mpf(delta)
    num = wh.stress_between_numeric(X - delta, X, mp.mpf(0), D, b, beta, abs(delta), **kw)
    xi = kt.doppler_xi(mp.mpf(beta))
    return num, kt.stress_between_kk(X - delta, X, mp.mpf(0), D, b, xi, num["n_terms"]).value


def _dev(num, S, name="eq4a", target=1):
    return max(abs(num[name][c].value / S - target) for c in COMPONENTS)


def _off_axis(rho_over_X, delta):
    b, X, delta = mp.mpf("0.01"), mp.mpf("0.25"), mp.mpf(delta)
    xi = kt.doppler_xi(mp.mpf("0.5"))
    rho = mp.mpf(rho_over_X) * X
    T = X - delta + (rho ** 2 / (2 * X) + rho ** 2 * xi / (2 * (D - X))) / (1 - xi)
    num = wh.stress_between_numeric(T, X, rho, D, b, "0.5", delta)
    return num, kt.stress_between_kk(T, X, rho, D, b, xi, num["n_terms"]).value


def test_origins():
    num, _ = _axis("0.5", "0.01", "0.25", "1e-3")
    assert num["eq4a"]["TT"].origin is Origin.COMPUTED
    xi = kt.doppler_xi(mp.mpf("0.5"))
    assert kt.stress_between_kk(mp.mpf(0), mp.mpf("0.25"), mp.mpf(0), D, mp.mpf("0.01"), xi, 3).origin is Origin.COMPUTED  # truncated series
    assert kt.stress_between_kk(mp.mpf(0), mp.mpf("0.25"), mp.mpf(0), D, mp.mpf("0.01"), xi, 3, [1]).origin is Origin.DERIVED  # one printed term


@pytest.mark.parametrize("beta", ["0.1", "0.5", "0.9"])
@pytest.mark.parametrize("b", ["0.01", "0.1"])
def test_eq4a_as_printed_reproduces_eq67_and_the_deviation_falls_with_delta(beta, b):
    near, S_near = _axis(beta, b, "0.25", "1e-6")
    far, S_far = _axis(beta, b, "0.25", "1e-4")
    assert _dev(near, S_near) < mp.mpf("0.41") * mp.mpf("1e-6")  # measured at most 0.404 delta
    assert abs(_dev(far, S_far) / _dev(near, S_near) / 100 - 1) < mp.mpf("1e-2")  # linear in delta


@pytest.mark.parametrize("beta", ["0.1", "0.5", "0.9"])
def test_operator_with_one_half_gives_half_of_eq67(beta):
    num, S = _axis(beta, "0.01", "0.25", "1e-6")
    assert _dev(num, S, "with_half", mp.mpf(1) / 2) < mp.mpf("0.21") * mp.mpf("1e-6")  # half of the above
    assert abs(num["eq4a"]["TT"].value / num["with_half"]["TT"].value - 2) < mp.mpf("1e-19")  # measured at most 2.0e-20


@pytest.mark.parametrize("beta", ["0.1", "0.5", "0.9"])
def test_midpoint_has_no_first_order_prefactor_term(beta):
    # d/dX of 1 / [X (D - X)] vanishes at X = D / 2, so T^TT and T^TX match eq. 67 down to the numerical error.
    num, S = _axis(beta, "0.01", "0.5", "1e-3")
    assert abs(num["eq4a"]["TT"].value / S - 1) < mp.mpf("1e-19")  # measured at most 4.0e-20
    assert abs(num["eq4a"]["TX"].value / S - 1) < mp.mpf("1e-19")


def test_finite_difference_error():
    a, _ = _axis("0.1", "0.01", "0.25", "1e-4")
    c, _ = _axis("0.1", "0.01", "0.25", "1e-4", split_ratio="1e-11", step_ratio="1e-12")
    assert abs(a["eq4a"]["TT"].value / c["eq4a"]["TT"].value - 1) < mp.mpf("1e-19")  # measured 4.0e-20


def test_series_truncation_is_below_the_numerical_error():
    a = wh.stress_between_numeric(mp.mpf("0.2499"), mp.mpf("0.25"), mp.mpf(0), D, mp.mpf("0.1"), "0.5", mp.mpf("1e-4"))
    c = wh.stress_between_numeric(mp.mpf("0.2499"), mp.mpf("0.25"), mp.mpf(0), D, mp.mpf("0.1"), "0.5", mp.mpf("1e-4"), tol="1e-40")
    assert c["n_terms"] > a["n_terms"]
    assert abs(a["eq4a"]["TT"].value / c["eq4a"]["TT"].value - 1) < mp.mpf("1e-32")  # measured 3.4e-33


def test_off_axis_corrections_scale_as_the_terms_the_article_ignores():
    big, S_big = _off_axis("1e-1", "1e-6")
    small, S_small = _off_axis("1e-2", "1e-6")
    assert abs(_dev(big, S_big) / _dev(small, S_small) / 100 - 1) < mp.mpf("0.05")  # (rho / X)^2, measured 0.00495 / 4.96e-5
    assert abs(big["eq4a"]["rest"].value / abs(S_big) / (small["eq4a"]["rest"].value / abs(S_small)) / 10 - 1) < mp.mpf("0.01")  # rho / X


def test_divergence_slope_and_sign_change_across_the_hypersurface():
    pts = [(mp.log(mp.mpf(d)), mp.log(abs(_axis("0.5", "0.01", "0.25", d)[0]["eq4a"]["TT"].value))) for d in ("1e-3", "1e-6")]
    assert abs((pts[1][1] - pts[0][1]) / (pts[1][0] - pts[0][0]) + 3) < mp.mpf("1e-5")  # measured 4.4e-6
    before, _ = _axis("0.5", "0.01", "0.25", "1e-3")
    after, _ = _axis("0.5", "0.01", "0.25", "-1e-3")
    assert before["eq4a"]["TT"].value < 0 < after["eq4a"]["TT"].value


# Dimensional analysis: every length scales with LAM.
LAM = mp.mpf(3)
TOL = mp.mpf("1e-75")
ARGS = (mp.mpf("0.2"), mp.mpf("0.25"), mp.mpf("0.01"), mp.mpf(0))
ARGS2 = (mp.mpf("0.19"), mp.mpf("0.26"), mp.mpf("0.02"), mp.mpf(0))
XI = kt.doppler_xi(mp.mpf("0.5"))
scaled = lambda v: [LAM * c for c in v]


def test_dimensions_of_eq61():
    ratio = kt.lambda_between(2, scaled(ARGS), scaled(ARGS2), LAM * D, XI) / kt.lambda_between(2, ARGS, ARGS2, D, XI)
    assert abs(ratio / LAM - 1) < TOL


def test_dimensions_of_eq66():
    b = mp.mpf("0.01")
    ratio = kt.hadamard_between(scaled(ARGS), scaled(ARGS2), LAM * D, LAM * b, XI, 3) / kt.hadamard_between(ARGS, ARGS2, D, b, XI, 3)
    assert abs(ratio * LAM ** 2 - 1) < TOL


def test_dimensions_of_eq67():
    T, X, rho, b = mp.mpf("0.2"), mp.mpf("0.25"), mp.mpf("0.01"), mp.mpf("0.01")
    a = kt.stress_between_kk(T, X, rho, D, b, XI, 3).value
    c = kt.stress_between_kk(LAM * T, LAM * X, LAM * rho, LAM * D, LAM * b, XI, 3).value
    assert abs(c / a * LAM ** 4 - 1) < TOL
