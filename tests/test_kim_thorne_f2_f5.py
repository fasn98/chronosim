"""F2 (defocusing) and F5 (polarized hypersurfaces) of Kim & Thorne.

Each tolerance sits a little above the measured value, quoted next to it.
"""

import mpmath as mp
import pytest

from chronosim import kim_thorne as kt
from chronosim import wormhole as wh
from chronosim.provenance import Origin

D, X = mp.mpf(1), mp.mpf("0.25")
BETAS = ("0.1", "0.5", "0.9")


@pytest.mark.parametrize("N", [1, 2, 3, 5])
def test_ray_transfer_reproduces_eq47_to_first_order_in_b_over_D(N):
    near = wh.ray_amplitude_throat(N, D, mp.mpf("1e-6")).value / kt.psi_throat(N, D, mp.mpf("1e-6")) - 1
    far = wh.ray_amplitude_throat(N, D, mp.mpf("1e-4")).value / kt.psi_throat(N, D, mp.mpf("1e-4")) - 1
    assert abs(near) <= (N - 1) * mp.mpf("1.01e-6")  # measured -(N - 1) b/D
    if N > 1:
        assert abs(far / near / 100 - 1) < mp.mpf("1e-3")


@pytest.mark.parametrize("N", [1, 2, 3, 5])
def test_ray_transfer_reproduces_eq64_to_first_order_in_b_over_D(N):
    xi = kt.doppler_xi(mp.mpf("0.5"))
    ratio = lambda b: wh.ray_amplitude_between(N, X, D, b).value / (kt.vleck_between(N, X, D, b, xi) / kt.zeta_between(N, X, D, xi)) - 1
    near, far = ratio(mp.mpf("1e-6")), ratio(mp.mpf("1e-4"))
    assert abs(near) < mp.mpf("7e-6")  # measured at most 6.7e-6
    assert abs(far / near / 100 - 1) < mp.mpf("1e-3")


def test_eq24_links_eq47_and_eq48():
    xi, b = kt.doppler_xi(mp.mpf("0.5")), mp.mpf("0.01")
    for N in (1, 2, 5):
        assert abs(kt.vleck_throat(N, D, b, xi) - kt.zeta(N, D, xi) * kt.psi_throat(N, D, b)) < mp.mpf("1e-75")


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("N", [1, 2, 5])
def test_eq57_is_the_exact_lapse_on_the_horizon(beta, N):
    xi = kt.doppler_xi(mp.mpf(beta))
    assert abs(kt.zeta_between_sum(N, X, D, xi) - kt.zeta_between(N, X, D, xi)) < mp.mpf("1e-78")
    assert abs(wh.exact_lapse_on_horizon(N, X, D, beta) / kt.zeta_between(N, X, D, xi) - 1) < mp.mpf("1e-79")  # measured 4.2e-81


def test_origins():
    assert wh.ray_amplitude_throat(2, D, mp.mpf("0.01")).origin is Origin.COMPUTED
    assert wh.root_HN_throat(1, mp.mpf("0.7"), mp.mpf("0.01"), "0.5").origin is Origin.COMPUTED


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("theta", ["0.2", "1.3"])
@pytest.mark.parametrize("N", [1, 2, 5, 20])
def test_roots_of_eq46_are_eq43(beta, theta, N):
    xi, b, theta = kt.doppler_xi(mp.mpf(beta)), mp.mpf("0.01"), mp.mpf(theta)
    assert abs(wh.root_HN_throat(N, theta, b, beta).value / kt.time_HN_throat(N, theta, b, xi) - 1) < mp.mpf("1e-78")  # measured 1.1e-80


@pytest.mark.parametrize("beta", BETAS)
@pytest.mark.parametrize("rho_over_X", ["1e-2", "1e-1"])
@pytest.mark.parametrize("N", [1, 2, 5, 20])
def test_roots_of_eq61_are_eq62(beta, rho_over_X, N):
    xi, rho = kt.doppler_xi(mp.mpf(beta)), mp.mpf(rho_over_X) * X
    T = kt.time_HN_between(N, X, rho, D, xi)
    assert abs((wh.root_HN_between(N, X, rho, D, beta).value - T) / (T - X)) < mp.mpf("1e-75")  # measured 1.8e-77


@pytest.mark.parametrize("beta", BETAS)
def test_hypersurfaces_approach_the_cauchy_horizon_as_xi_to_the_N(beta):
    xi, b, theta = kt.doppler_xi(mp.mpf(beta)), mp.mpf("0.01"), mp.mpf("0.7")
    for N in (10, 30):
        gap = wh.root_HN_throat(N, theta, b, beta).value / kt.time_Hminus_throat(theta, b) - 1
        assert abs(gap / (2 * xi ** N / (1 - xi ** N)) - 1) < mp.mpf("1e-60")


def test_hypersurfaces_are_nested_in_decreasing_N():
    xi, b, theta = kt.doppler_xi(mp.mpf("0.5")), mp.mpf("0.01"), mp.mpf("0.7")
    times = [kt.time_HN_throat(N, theta, b, xi) for N in range(1, 8)]
    assert all(times[i] > times[i + 1] > kt.time_Hminus_throat(theta, b) for i in range(6))


# Dimensional analysis: every length scales with LAM.
LAM = mp.mpf(3)


def test_dimensions():
    xi, b, theta, rho = kt.doppler_xi(mp.mpf("0.5")), mp.mpf("0.01"), mp.mpf("0.7"), mp.mpf("0.01")
    assert abs(kt.psi_throat(3, LAM * D, LAM * b) * LAM / kt.psi_throat(3, D, b) - 1) < mp.mpf("1e-75")  # 1 / L
    assert abs(kt.vleck_throat(3, LAM * D, LAM * b, xi) / kt.vleck_throat(3, D, b, xi) - 1) < mp.mpf("1e-75")  # 1
    assert abs(kt.vleck_between(3, LAM * X, LAM * D, LAM * b, xi) / kt.vleck_between(3, X, D, b, xi) - 1) < mp.mpf("1e-75")
    assert abs(kt.zeta_between(3, LAM * X, LAM * D, xi) / kt.zeta_between(3, X, D, xi) / LAM - 1) < mp.mpf("1e-75")
    assert abs(kt.time_HN_throat(3, theta, LAM * b, xi) / kt.time_HN_throat(3, theta, b, xi) / LAM - 1) < mp.mpf("1e-75")
    assert abs(kt.time_Hminus_throat(theta, LAM * b) / kt.time_Hminus_throat(theta, b) / LAM - 1) < mp.mpf("1e-75")
    assert abs(kt.time_HN_between(3, LAM * X, LAM * rho, LAM * D, xi) / kt.time_HN_between(3, X, rho, D, xi) / LAM - 1) < mp.mpf("1e-75")
    assert abs(kt.time_Hminus_between(LAM * X, LAM * rho) / kt.time_Hminus_between(X, rho) / LAM - 1) < mp.mpf("1e-75")
