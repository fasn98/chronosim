"""F4: Kim & Thorne eqs. 46 and 50 -> eq. 51, without using eq. 51 in the numerics.

Each tolerance sits a little above the measured value, quoted next to it.
"""

import mpmath as mp
import pytest

from chronosim import kim_thorne as kt
from chronosim import wormhole as wh
from chronosim.provenance import Origin

D = mp.mpf(1)


def _point(beta, b, theta, delta_over_b, phi="0"):
    b, theta, phi = mp.mpf(b), mp.mpf(theta), mp.mpf(phi)
    xi = kt.doppler_xi(mp.mpf(beta))
    delta = mp.mpf(delta_over_b) * b
    return b, theta, phi, xi, delta, kt.hypersurface_time_throat(1, theta, b, xi) - delta


def _dev(beta, b, theta, delta_over_b, variant, name="eq4a", phi="0"):
    b, theta, phi, xi, delta, T = _point(beta, b, theta, delta_over_b, phi)
    num = wh.stress_throat_numeric(kt.throat_point(T, theta, phi, b), D, b, beta, delta)
    E = kt.stress_throat(T, theta, phi, D, b, xi, num["n_terms"], variant)
    target = mp.mpf(1) if name == "eq4a" else mp.mpf(1) / 2
    scale = max(abs(E[m][c].value) for m in range(4) for c in range(4))
    return max(abs(num[name][m][c].value - target * E[m][c].value) for m in range(4) for c in range(4)) / (target * scale)


def test_origins():
    b, theta, phi, xi, delta, T = _point("0.5", "0.01", "0.7", "1e-3")
    assert kt.stress_throat(T, theta, phi, D, b, xi, 3, "transverse")[0][0].origin is Origin.COMPUTED  # truncated series
    assert kt.stress_throat(T, theta, phi, D, b, xi, 3, "transverse", [2])[0][0].origin is Origin.DERIVED  # one printed term
    num = wh.stress_throat_numeric(kt.throat_point(T, theta, phi, b), D, b, "0.5", delta)
    assert num["eq4a"][0][0].origin is Origin.COMPUTED


@pytest.mark.parametrize("beta", ["0.1", "0.5", "0.9"])
@pytest.mark.parametrize("theta", ["0.2", "0.7", "1.3"])
def test_eq4a_as_printed_reproduces_eq51_with_the_transverse_rho_hat(beta, theta):
    assert _dev(beta, "0.01", theta, "1e-6", "transverse") < mp.mpf("8e-7")  # measured at most 7.8e-7
    ratio = _dev(beta, "0.01", theta, "1e-4", "transverse") / _dev(beta, "0.01", theta, "1e-6", "transverse")
    assert abs(ratio / 100 - 1) < mp.mpf("1e-2")  # linear in delta


def test_operator_with_one_half_gives_half_of_eq51():
    assert _dev("0.5", "0.01", "0.7", "1e-6", "transverse", "with_half") < mp.mpf("3e-7")  # measured 2.6e-7


def test_also_off_the_y_axis():
    assert _dev("0.5", "0.01", "0.7", "1e-6", "transverse", phi="0.6") < mp.mpf("3e-7")  # measured 2.3e-7


@pytest.mark.parametrize("theta,expected", [("0.7", "0.09"), ("1.3", "0.57")])
def test_literal_denominator_of_eq52_does_not_reproduce_eq51(theta, expected):
    # Its deviation does not fall with delta. Measured 0.0916 and 0.5779 at delta / b = 1e-6.
    near = _dev("0.5", "0.01", theta, "1e-6", "literal")
    far = _dev("0.5", "0.01", theta, "1e-4", "literal")
    assert near > mp.mpf(expected)
    assert abs(near / far - 1) < mp.mpf("1e-3")


def test_literal_rho_hat_is_not_unit_on_the_throat():
    b, theta = mp.mpf("0.01"), mp.mpf("0.7")
    _, X, Y, Z = kt.throat_point(mp.mpf(0), theta, mp.mpf(0), b)
    assert abs(kt.rho_hat(X, Y, Z, "literal")[2] - mp.cos(theta / 2)) < mp.mpf("1e-75")
    assert abs(kt.rho_hat(X, Y, Z, "transverse")[2] - 1) < mp.mpf("1e-75")


@pytest.mark.parametrize("beta,theta,phi", [("0.5", "0.7", "0"), ("0.9", "1.3", "0"), ("0.5", "0.7", "0.6")])
def test_residual_is_the_second_derivative_term_that_eq51_omits(beta, theta, phi):
    # With rho = (Y, Z), eq. 4b has -1/3 d_i d_j of 1 / lambda_N(x, x'), which contains
    # +1/3 (d_i d_j lambda_N) / lambda_N^2. Eq. 51 keeps only the 1 / lambda^3 terms.
    b, theta, phi, xi, delta, T = _point(beta, "0.01", theta, "1e-4", phi)
    P = kt.throat_point(T, theta, phi, b)
    num = wh.stress_throat_numeric(P, D, b, beta, delta)
    n = num["n_terms"]
    E = kt.stress_throat(T, theta, phi, D, b, xi, n, "transverse")
    transverse = (P[2], P[3])
    root = mp.sqrt(b ** 2 - transverse[0] ** 2 - transverse[1] ** 2)
    scale = max(abs(E[m][c].value) for m in range(4) for c in range(4))
    for i in range(2):
        for j in range(2):
            hessian = 2 * ((1 if i == j else 0) / root + transverse[i] * transverse[j] / root ** 3)
            predicted = sum(
                xi / (4 * mp.pi ** 2 * D) * (b * xi / (2 * D)) ** (N - 1) * hessian / 3 / kt.lambda_throat(N, P, P, b, xi) ** 2
                for N in range(1, n + 1)
            )
            residual = num["eq4a"][2 + i][2 + j].value - E[2 + i][2 + j].value
            assert abs(residual - predicted) / scale < mp.mpf("1e-8")
    for m in range(4):
        for c in range(4):
            if m < 2 or c < 2:
                assert abs(num["eq4a"][m][c].value - E[m][c].value) / scale < mp.mpf("1e-17")  # measured at most 5.1e-19


def test_divergence_slope_and_sign_change_across_the_hypersurface():
    def tt(delta_over_b):
        b, theta, phi, xi, delta, T = _point("0.5", "0.01", "0.7", delta_over_b)
        return wh.stress_throat_numeric(kt.throat_point(T, theta, phi, b), D, b, "0.5", abs(delta))["eq4a"][0][0].value

    near, far = tt("1e-6"), tt("1e-3")
    assert abs(mp.log(abs(near / far)) / mp.log(mp.mpf("1e-3")) + 3) < mp.mpf("1e-11")  # measured 9.6e-13
    assert tt("1e-4") < 0 < tt("-1e-4")


# Dimensional analysis: every length scales with LAM.
LAM = mp.mpf(3)
TOL = mp.mpf("1e-72")
B = mp.mpf("0.01")
XI = kt.doppler_xi(mp.mpf("0.5"))
A = [mp.mpf("-0.001"), mp.mpf("-0.0002"), mp.mpf("0.005"), mp.mpf("0.002")]
C = [mp.mpf("-0.0011"), mp.mpf("-0.0002"), mp.mpf("0.005"), mp.mpf("0.002")]
scaled = lambda v: [LAM * c for c in v]


def test_dimensions_of_eq46():
    ratio = kt.lambda_throat(2, scaled(A), scaled(C), LAM * B, XI) / kt.lambda_throat(2, A, C, B, XI)
    assert abs(ratio / LAM - 1) < TOL


def test_dimensions_of_eq50():
    ratio = kt.hadamard_throat(scaled(A), scaled(C), LAM * D, LAM * B, XI, 3) / kt.hadamard_throat(A, C, D, B, XI, 3)
    assert abs(ratio * LAM ** 2 - 1) < TOL


@pytest.mark.parametrize("variant", ["literal", "transverse"])
def test_dimensions_of_eq51_and_eq52(variant):
    theta, phi, T = mp.mpf("0.7"), mp.mpf("0.6"), mp.mpf("-0.001")
    a = kt.stress_throat(T, theta, phi, D, B, XI, 3, variant)
    c = kt.stress_throat(LAM * T, theta, phi, LAM * D, LAM * B, XI, 3, variant)
    assert abs(c[0][0].value / a[0][0].value * LAM ** 4 - 1) < TOL
    assert abs(c[2][2].value / a[2][2].value * LAM ** 4 - 1) < TOL


def test_eq52_chain_rule_identity():
    # rho_hat = d/drho with rho = sqrt(Y^2 + Z^2) (eq. 39): d/drho = (Y d/dY + Z d/dZ) / sqrt(Y^2 + Z^2).
    # Checked on f(Y, Z) = g(rho) with g = sin. The printed denominator sqrt(X^2 + Y^2) fails it.
    mp.mp.dps = 50
    X, Y, Z = mp.mpf("-1"), mp.mpf("0.4"), mp.mpf("0.25")
    rho = mp.sqrt(Y ** 2 + Z ** 2)
    f = lambda y, z: mp.sin(mp.sqrt(y ** 2 + z ** 2))
    numerator = Y * mp.diff(lambda y: f(y, Z), Y) + Z * mp.diff(lambda z: f(Y, z), Z)
    assert abs(numerator / mp.sqrt(Y ** 2 + Z ** 2) - mp.cos(rho)) < mp.mpf("1e-40")
    assert abs(numerator / mp.sqrt(X ** 2 + Y ** 2) - mp.cos(rho)) > mp.mpf("0.1")
