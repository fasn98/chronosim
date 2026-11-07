"""F1: the covering space of Kim & Thorne, built independently from their eqs. 25 and 28.

Each tolerance sits a little above the measured error, quoted next to it. That error
is rounding at dps = 80, amplified by xi^-N where it applies.
"""

import mpmath as mp
import pytest

from chronosim import kim_thorne as kt
from chronosim import wormhole as wh
from chronosim.provenance import Origin

D = mp.mpf(1)
BETAS = ("0.1", "0.5", "0.9")
NS = (1, 2, 3, 5, 10)
EPS = ("1e-2", "1e-4", "1e-6", "1e-8")


def _pair(eps):
    eps = mp.mpf(eps)
    return (-eps, 2 * eps), (-2 * eps, eps + eps / 3)


@pytest.mark.parametrize("beta", BETAS)
def test_map_fixes_the_covering_space_origin(beta):
    beta = mp.mpf(beta)
    fixed = D / (1 - kt.doppler_xi(beta))
    T, X = wh.throat_map((fixed, fixed), D, beta)
    assert max(abs(T - fixed), abs(X - fixed)) < mp.mpf("1e-79")  # measured 2.1e-81


@pytest.mark.parametrize("beta", BETAS)
def test_map_sends_mouth_1_to_mouth_2_at_equal_proper_time(beta):
    # Eqs. 25b, 26, 28: the image of (tau, 0) has X' = 0 and T' = tau.
    beta, tau = mp.mpf(beta), mp.mpf("0.7")
    gamma = 1 / mp.sqrt(1 - beta ** 2)
    T, X = wh.throat_map((tau, mp.mpf(0)), D, beta)
    assert abs(T - (D + gamma * tau)) < mp.mpf("1e-79")
    assert abs(X - (D - beta * gamma * tau)) < mp.mpf("1e-79")
    assert abs(gamma * (X - D) + beta * gamma * (T - D)) < mp.mpf("1e-79")  # measured 1.3e-82
    assert abs(gamma * (T - D) + beta * gamma * (X - D) - tau) < mp.mpf("1e-79")  # measured 2.1e-81


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("n", NS)
def test_eq33_copies(beta, n):
    beta = mp.mpf(beta)
    xi = kt.doppler_xi(beta)
    event = (mp.mpf("0.123"), mp.mpf("0.456"))
    U0, V0 = kt.null_coords(kt.to_covering(event, D, xi))
    U, V = kt.null_coords(kt.to_covering(wh.copy_event(event, n, D, beta), D, xi))
    Un, Vn = kt.copy_n(U0, V0, n, xi)
    assert abs(U / Un - 1) < mp.mpf("1e-68")  # measured at most 3.9e-70
    assert abs(V / Vn - 1) < mp.mpf("1e-68")


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("n", (1, 2, 3))
@pytest.mark.parametrize("tau", ("-1", "0.7"))
def test_eq31_as_printed_differs_from_the_map_by_a_sign(beta, n, tau):
    # The copy T_n of the throat, from eqs. 25, 28 and 29, satisfies eq. 31 only
    # with the opposite sign of its coefficient. Both transcriptions read "+".
    beta = mp.mpf(beta)
    xi = kt.doppler_xi(beta)
    T, X = kt.to_covering(wh.copy_event((mp.mpf(tau), mp.mpf(0)), n, D, beta), D, xi)
    c = xi ** n * D / (1 - xi)
    coefficient = (1 + xi ** (2 * n)) / (1 - xi ** (2 * n))
    assert abs((T + c) / (X + c) / coefficient + 1) < mp.mpf("1e-65")
    assert abs(kt.throat_copy_residual((T, X), n, D, xi)) > mp.mpf("1e-3")


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("n", (1, 2, 3))
def test_throat_copies_have_negative_slope_as_in_fig_5(beta, n):
    # Fig. 5 draws T_1, T_2, ... with negative slope in the X-T plane. Eq. 31 as printed
    # has a positive slope for 0 < xi < 1, since its coefficient is positive.
    beta = mp.mpf(beta)
    xi = kt.doppler_xi(beta)
    Ta, Xa = kt.to_covering(wh.copy_event((mp.mpf(-1), mp.mpf(0)), n, D, beta), D, xi)
    Tb, Xb = kt.to_covering(wh.copy_event((mp.mpf(1), mp.mpf(0)), n, D, beta), D, xi)
    assert (Tb - Ta) / (Xb - Xa) < 0
    assert (1 + xi ** (2 * n)) / (1 - xi ** (2 * n)) > 0


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("tau", ("-1", "0", "0.7"))
def test_copy_zero_of_the_throat_is_at_the_declared_position(beta, tau):
    # Text before eq. 31: T_0 is at X~ = -D / (1 - xi).
    xi = kt.doppler_xi(mp.mpf(beta))
    _, X = kt.to_covering((mp.mpf(tau), mp.mpf(0)), D, xi)
    assert abs(X + D / (1 - xi)) < mp.mpf("1e-78")


@pytest.mark.parametrize("beta", BETAS)
def test_horizon_crossings_reproduce_the_spacings_of_fig_5(beta):
    # tau = 0 lies on the horizon. Copy n crosses it at X~ = T~ = -xi^n D / (1 - xi),
    # so neighbouring crossings are D, D xi, D xi^2 apart.
    beta = mp.mpf(beta)
    xi = kt.doppler_xi(beta)
    crossing = [kt.to_covering(wh.copy_event((mp.mpf(0), mp.mpf(0)), n, D, beta), D, xi) for n in range(5)]
    for n, (T, X) in enumerate(crossing):
        assert abs(T - X) < mp.mpf("1e-77")
        assert abs(X + xi ** n * D / (1 - xi)) < mp.mpf("1e-77")
    for n in range(4):
        assert abs(crossing[n + 1][1] - crossing[n][1] - D * xi ** n) < mp.mpf("1e-77")


@pytest.mark.parametrize("sign", (1, -1))
@pytest.mark.parametrize("n", (1, 2, 3))
def test_position_and_crossings_do_not_distinguish_the_sign_of_eq31(sign, n):
    # Written as (1 - xi^2n)(T~ + c) = sign (1 + xi^2n)(X~ + c), so that n = 0 is defined. Both
    # signs give X~ = -c at n = 0 and the same horizon crossing. Only the slope tells them apart.
    xi = kt.doppler_xi(mp.mpf("0.5"))
    residual = lambda T, X, m: (
        (1 - xi ** (2 * m)) * (T + xi ** m * D / (1 - xi))
        - sign * (1 + xi ** (2 * m)) * (X + xi ** m * D / (1 - xi))
    )
    assert abs(residual(mp.mpf("0.3"), -D / (1 - xi), 0)) < mp.mpf("1e-78")
    crossing = -xi ** n * D / (1 - xi)
    assert abs(residual(crossing, crossing, n)) < mp.mpf("1e-78")


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("N", NS)
def test_eq35_differs_from_the_exact_lapse_by_first_order_terms(beta, N):
    # V_N(x) - V_0(x') = zeta_N + xi^N (X + T) / 2 - (X' + T') / 2, from eqs. 29, 32 and 33.
    # The article says only "very nearly the same"; this is the source of the linear error of eq. 38.
    beta = mp.mpf(beta)
    xi = kt.doppler_xi(beta)
    x, xp = _pair("1e-4")
    _, V = kt.null_coords(kt.to_covering(wh.copy_event(x, N, D, beta), D, xi))
    _, Vp = kt.null_coords(kt.to_covering(xp, D, xi))
    expected = kt.zeta(N, D, xi) + xi ** N * (x[1] + x[0]) / 2 - (xp[1] + xp[0]) / 2
    assert abs(V - Vp - expected) < mp.mpf("1e-76")


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("N", NS)
def test_eq36_closed_form_is_the_sum(beta, N):
    xi = kt.doppler_xi(mp.mpf(beta))
    assert abs(kt.zeta(N, D, xi) - sum(D * xi ** n for n in range(N))) < mp.mpf("1e-78")


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("N", NS)
@pytest.mark.parametrize("eps", EPS)
def test_eqs_34_and_37_are_exact(beta, N, eps):
    beta = mp.mpf(beta)
    xi = kt.doppler_xi(beta)
    x, xp = _pair(eps)
    U, V = kt.null_coords(kt.to_covering(wh.copy_event(x, N, D, beta), D, xi))
    Up, Vp = kt.null_coords(kt.to_covering(xp, D, xi))
    exact = wh.exact_interval(N, x, xp, D, beta)
    assert abs((V - Vp) * (U - Up) / exact - 1) < mp.mpf("1e-68")  # measured at most 1.5e-70
    assert abs((U - Up) / kt.delta_U(N, x, xp, xi) - 1) < mp.mpf("1e-68")  # measured at most 1.5e-70


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("N", NS)
def test_eq38_error_falls_linearly_with_the_distance_to_the_throat(beta, N):
    beta = mp.mpf(beta)
    err = [abs(kt.sigma_N(N, *_pair(e), D, beta).value / wh.exact_interval(N, *_pair(e), D, beta) - 1) for e in EPS]
    assert max(e / mp.mpf(s) for e, s in zip(err, EPS)) < mp.mpf("0.8")  # measured at most 0.786
    assert abs(err[2] / err[3] / 100 - 1) < mp.mpf("1e-5")


def test_origins():
    x, xp = _pair("1e-4")
    assert kt.sigma_N(1, x, xp, D, "0.5").origin is Origin.DERIVED
    assert wh.exact_sigma_N(1, x, xp, D, "0.5").origin is Origin.COMPUTED


# Dimensional analysis: with c = G = hbar = 1, every length scales with lam, so each
# formula must scale with the power of lam that its dimension requires.
LAM = mp.mpf(3)
TOL = mp.mpf("1e-75")


def _scaled_ratio(f, *lengths, power):
    return abs(f(*[LAM * v for v in lengths]) / f(*lengths) / LAM ** power - 1)


def test_dimensions_of_the_flat_hadamard_function():  # eqs. 3, 17: 1 / L^2 with sigma ~ L^2
    assert abs(kt.hadamard_flat(LAM ** 2 * mp.mpf(2)) / kt.hadamard_flat(mp.mpf(2)) * LAM ** 2 - 1) < TOL


def test_eq30_is_dimensionless_and_below_one():
    assert 0 < kt.doppler_xi(mp.mpf("0.5")) < 1


def test_dimensions_of_eq29():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    assert _scaled_ratio(lambda T, X, d: kt.to_covering((T, X), d, xi)[0], mp.mpf("0.1"), mp.mpf("0.2"), D, power=1) < TOL


def test_dimensions_of_eq32_and_eq33():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    U = lambda T, X: kt.copy_n(*kt.null_coords((T, X)), 2, xi)[0]
    V = lambda T, X: kt.copy_n(*kt.null_coords((T, X)), 2, xi)[1]
    assert _scaled_ratio(U, mp.mpf("0.1"), mp.mpf("0.2"), power=1) < TOL
    assert _scaled_ratio(V, mp.mpf("0.1"), mp.mpf("0.2"), power=1) < TOL


def test_dimensions_of_eq31():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    f = lambda T, X, d: kt.throat_copy_residual((T, X), 1, d, xi)
    assert _scaled_ratio(f, mp.mpf("0.1"), mp.mpf("0.2"), D, power=1) < TOL


def test_dimensions_of_eq36_and_eq37():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    assert _scaled_ratio(lambda d: kt.zeta(3, d, xi), D, power=1) < TOL
    f = lambda T, X, Tp, Xp: kt.delta_U(3, (T, X), (Tp, Xp), xi)
    assert _scaled_ratio(f, mp.mpf("-0.1"), mp.mpf("0.2"), mp.mpf("-0.2"), mp.mpf("0.1"), power=1) < TOL


def test_dimensions_of_eq38_and_the_exact_interval():
    args = (mp.mpf("-0.1"), mp.mpf("0.2"), mp.mpf("-0.2"), mp.mpf("0.1"), D)
    beta = mp.mpf("0.5")
    f = lambda T, X, Tp, Xp, d: kt.sigma_N(3, (T, X), (Tp, Xp), d, beta).value
    g = lambda T, X, Tp, Xp, d: wh.exact_interval(3, (T, X), (Tp, Xp), d, beta)
    assert _scaled_ratio(f, *args, power=2) < TOL
    assert _scaled_ratio(g, *args, power=2) < TOL


def test_dimensions_of_eq4_on_the_rindler_vacuum():
    from chronosim.misner import hadamard_image

    F = wh.regularized_with_eq3(lambda P, Q: hadamard_image(0, P, Q, 1))
    rho = lambda xi: wh.eq4_stress(F, [0, xi, 0, 0], delta=xi * mp.mpf("1e-10"), h=xi * mp.mpf("1e-11"))[0][0]
    assert _scaled_ratio(rho, mp.mpf(1), power=-4) < mp.mpf("1e-15")
