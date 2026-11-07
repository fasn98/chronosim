"""F6: four-momentum by quadrature of eqs. 51 and 67, against eqs. 53, 54, 68 and 69.

Planck units, l_P = m_P = 1. The quadratures run at kt_quad_dps = 15, so agreement
is limited to about 1e-8 where the comparison is exact. Each tolerance sits above the
measured value, quoted next to it. These tests take about two minutes.
"""

import mpmath as mp
import pytest

from chronosim import kim_thorne as kt
from chronosim import wormhole as wh

D, B, X = mp.mpf(1000), mp.mpf(10), mp.mpf(250)


@pytest.mark.parametrize("beta", ["0.1", "0.5", "0.9"])
@pytest.mark.parametrize("N", [1, 3])
@pytest.mark.parametrize("rho_over_X", ["0", "1e-2"])
def test_eq69_and_its_independence_of_rho(beta, N, rho_over_X):
    xi = kt.doppler_xi(mp.mpf(beta))
    numeric = wh.momentum_per_area_numeric(N, X, mp.mpf(rho_over_X) * X, D, B, beta).value
    assert abs(numeric / kt.momentum_per_area_N(N, X, D, B, xi).value - 1) < mp.mpf("1e-7")  # measured at most 2.0e-8


def test_eq68_for_the_first_order():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    numeric = wh.momentum_plane_numeric(mp.mpf(1), X, D, B, "0.5", orders=[1]).value
    printed = -(B / (2 * D)) * (1 + 4 * xi + xi ** 2) / (1 - xi) ** 2 / (6 * mp.pi)
    assert abs(numeric / printed - 1) < mp.mpf("1e-9")  # measured 5.5e-11


def test_eq54_and_eq68_are_the_same_expression():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    assert kt.momentum_horizon(mp.mpf(1), D, B, xi, 5).value < 0


def test_eq54_is_the_leading_term_as_delta_T_over_b_falls():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    single = -(B / (2 * D)) * (1 + 4 * xi + xi ** 2) / (1 - xi) ** 2 / (6 * mp.pi)
    dev = [abs(wh.momentum_throat_numeric(mp.mpf(dt), D, B, "0.5", lambda th: kt.time_Hminus_throat(th, B), orders=[1]).value
               / (single / mp.mpf(dt)) - 1) for dt in ("0.1", "0.01")]
    assert dev[1] < mp.mpf("0.01")  # measured 0.0091
    assert dev[1] < dev[0] / 5  # measured 0.079 -> 0.0091


def test_eq53_coefficient_of_the_logarithm():
    xi = kt.doppler_xi(mp.mpf("0.5"))
    upper = lambda th: kt.time_HN_throat(1, th, B, xi)
    v4 = wh.momentum_throat_numeric(mp.mpf(1), D, B, "0.5", upper, orders=[1], thickness="1e-4").value
    v6 = wh.momentum_throat_numeric(mp.mpf(1), D, B, "0.5", upper, orders=[1], thickness="1e-6").value
    coefficient = kt.momentum_throat_N(1, mp.mpf(1), B / mp.e, D, B, xi).value
    assert abs((v6 - v4) / (2 * mp.log(10)) / coefficient - 1) < mp.mpf("1e-7")  # measured 1.2e-9


LAM = mp.mpf(3)


def test_dimensions():
    # P ~ 1 / length (m_P l_P / Delta T) and dP/dA ~ 1 / length^2 with hbar = 1.
    xi = kt.doppler_xi(mp.mpf("0.5"))
    assert abs(kt.momentum_horizon(LAM, LAM * D, LAM * B, xi, 4).value * LAM / kt.momentum_horizon(1, D, B, xi, 4).value - 1) < mp.mpf("1e-75")
    a = kt.momentum_throat_N(2, mp.mpf(1), mp.mpf("0.001"), D, B, xi).value
    c = kt.momentum_throat_N(2, LAM, LAM * mp.mpf("0.001"), LAM * D, LAM * B, xi).value
    assert abs(c * LAM / a - 1) < mp.mpf("1e-75")
    assert abs(kt.momentum_per_area_N(2, LAM * X, LAM * D, LAM * B, xi).value * LAM / kt.momentum_per_area_N(2, X, D, B, xi).value - 1) < mp.mpf("1e-75")
