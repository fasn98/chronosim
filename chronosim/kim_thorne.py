"""Formulas transcribed from Kim & Thorne (1991). See docs/equation_transcription.md.

Conventions (their Sec. I): c = G = hbar = 1, signs of Misner, Thorne & Wheeler.
Their xi is the inverse Doppler blueshift of eq. 30, not the Rindler coordinate
of ``misner.py``. Events in the T-X plane are (T, X).
"""

from __future__ import annotations

import mpmath as mp

from .provenance import Origin, Quantity, heuristic
from .sources import KIM_THORNE_1991


def hadamard_flat(sigma0):
    """Flat-space Hadamard function 1 / (4 pi^2 sigma), eqs. 3 and 17.

    sigma is half the square of the separation, the N = 0 term of eq. 20.
    """
    return 1 / (4 * mp.pi ** 2 * sigma0)


def doppler_xi(beta):
    """xi = [(1 - beta) / (1 + beta)]^(1/2), eq. 30."""
    return mp.sqrt((1 - beta) / (1 + beta))


def to_covering(event, D, xi):
    """(T, X) in copy 0 -> covering-space (T~, X~), eq. 29."""
    T, X = event
    shift = D / (1 - xi)
    return T - shift, X - shift


def null_coords(event_tilde):
    """(T~, X~) -> (U, V) = (X~ - T~, (X~ + T~) / 2), eq. 32."""
    T, X = event_tilde
    return X - T, (X + T) / 2


def copy_n(U0, V0, n, xi):
    """Copy n of an event at (U0, V0) in copy 0: U_n = xi^-n U0, V_n = xi^n V0, eq. 33."""
    return xi ** -n * U0, xi ** n * V0


def throat_copy_residual(event_tilde, n, D, xi):
    """Left minus right side of eq. 31 as printed. Zero on copy T_n if eq. 31 holds."""
    T, X = event_tilde
    c = xi ** n * D / (1 - xi)
    return (T + c) - (1 + xi ** (2 * n)) / (1 - xi ** (2 * n)) * (X + c)


def zeta(N, D, xi):
    """Total affine parameter lapse zeta_N = D (1 - xi^N) / (1 - xi), eq. 36."""
    return D * (1 - xi ** N) / (1 - xi)


def delta_U(N, x, xp, xi):
    """U_N(x) - U_0(x') = (X - T) xi^-N - (X' - T'), eq. 37."""
    T, X = x
    Tp, Xp = xp
    return (X - T) * xi ** -N - (Xp - Tp)


def sigma_N(N, x, xp, D, beta) -> Quantity:
    """Nth geodetic interval on the symmetry axis, sigma_N = zeta_N [(X - T) xi^-N - (X' - T')], eq. 38."""
    xi = doppler_xi(mp.mpf(beta))
    return Quantity(zeta(N, D, xi) * delta_U(N, x, xp, xi), Origin.DERIVED, KIM_THORNE_1991.ref("eq. 38"))


def _series_origin(orders, n_terms):
    """DERIVED for a single printed term; COMPUTED for a sum truncated at kt_series_tol (an image sum)."""
    count = n_terms if orders is None else len(list(orders))
    return Origin.DERIVED if count == 1 else Origin.COMPUTED


SERIES_TOL = heuristic("kt_series_tol", mp.mpf("1e-30"), "the eqs. 66-67 series stops when (b xi / 2D)^N is below this")


def series_terms(D, b, xi, tol=None):
    """Number of terms N with (b xi / 2D)^N above the tolerance."""
    tol = SERIES_TOL.value if tol is None else mp.mpf(tol)
    return int(mp.ceil(mp.log(tol) / mp.log(b * xi / (2 * D))))


def lambda_between(N, x, xp, D, xi):
    """lambda_N(x, x') between the mouths, eq. 61. Events are (T, X, Y, Z)."""
    T, X, Y, Z = x
    Tp, Xp, Yp, Zp = xp
    return (
        (X - T) - (Xp - Tp) * xi ** N
        + (Y ** 2 + Z ** 2) / (2 * X)
        + (Yp ** 2 + Zp ** 2) / (2 * (D - Xp)) * xi ** N
    )


def hadamard_between(x, xp, D, b, xi, n_terms):
    """Regularized Hadamard function between the mouths, eq. 66.

    Read literally: the prefactor uses the X of the first argument, not a symmetrized one.
    """
    X = x[1]
    prefactor = D / (4 * mp.pi ** 2 * X * (D - X))
    return sum(
        prefactor * (b * xi / (2 * D)) ** N
        * (1 / lambda_between(N, x, xp, D, xi) + 1 / lambda_between(N, xp, x, D, xi))
        for N in range(1, n_terms + 1)
    )


def stress_between_kk(T, X, rho, D, b, xi, n_terms, orders=None) -> Quantity:
    """Coefficient of k^mu k^nu in T^mu^nu between the mouths, eq. 67, as printed.

    ``orders`` restricts the sum to the given N (default 1..n_terms).
    """
    total = mp.mpf(0)
    for N in (range(1, n_terms + 1) if orders is None else orders):
        a = xi ** N
        denominator = (X - T) * (1 - a) + rho ** 2 / (2 * X) + rho ** 2 / (2 * (D - X)) * a
        total += (b * xi / (2 * D)) ** N * (1 + 4 * a + a ** 2) / denominator ** 3
    origin = _series_origin(orders, n_terms)
    note = "" if origin is Origin.DERIVED else ", series truncated at kt_series_tol"
    return Quantity(-D / (6 * mp.pi ** 2 * X * (D - X)) * total, origin, KIM_THORNE_1991.ref("eq. 67") + note)


def lambda_throat(N, x, xp, b, xi):
    """lambda_N(x, x') near the throat, eq. 46. Events are (T, X, Y, Z); rho follows eq. 39."""
    T, X, Y, Z = x
    Tp, Xp = xp[0], xp[1]
    rho2 = Y ** 2 + Z ** 2
    return 2 * (b - mp.sqrt(b ** 2 - rho2)) + X - T - (Xp - Tp) * xi ** N


def hadamard_throat(x, xp, D, b, xi, n_terms):
    """Regularized Hadamard function near the throat, eq. 50."""
    prefactor = xi / (4 * mp.pi ** 2 * D)
    return sum(
        prefactor * (b * xi / (2 * D)) ** (N - 1)
        * (1 / lambda_throat(N, x, xp, b, xi) + 1 / lambda_throat(N, xp, x, b, xi))
        for N in range(1, n_terms + 1)
    )


RHO_HAT_DENOMINATORS = {
    "literal": lambda X, Y, Z: mp.sqrt(X ** 2 + Y ** 2),  # eq. 52 as printed (checked by a human reader)
    "transverse": lambda X, Y, Z: mp.sqrt(Y ** 2 + Z ** 2),  # the rho of eq. 39
}


def rho_hat(X, Y, Z, variant):
    """(rho^T, rho^X, rho^Y, rho^Z) from eq. 52, with the chosen denominator."""
    n = RHO_HAT_DENOMINATORS[variant](X, Y, Z)
    return [mp.mpf(0), mp.mpf(0), Y / n, Z / n]


def throat_point(T, theta, phi, b):
    """(T, X, Y, Z) on the inner face of the throat, eq. 42: rho = b sin theta, X + b = b cos theta."""
    rho = b * mp.sin(theta)
    return [T, b * mp.cos(theta) - b, rho * mp.cos(phi), rho * mp.sin(phi)]


def hypersurface_time_throat(N, theta, b, xi):
    """T where the denominator of eq. 51 vanishes, for the given N (the same as eq. 43)."""
    a = xi ** N
    return (1 + a) / (1 - a) * b * (1 - mp.cos(theta))


def stress_throat(T, theta, phi, D, b, xi, n_terms, variant, orders=None, rows=range(4)) -> list:
    """T^mu^nu on the throat, eq. 51, as a 4x4 list of DERIVED quantities, with k from eq. 52.

    ``orders`` restricts the sum to the given N (default 1..n_terms); ``rows`` limits the
    rows computed (the others stay zero).
    """
    point = throat_point(T, theta, phi, b)
    k = [mp.mpf(1), mp.mpf(1), mp.mpf(0), mp.mpf(0)]
    rh = rho_hat(point[1], point[2], point[3], variant)
    tan = mp.tan(theta)
    out = [[mp.mpf(0)] * 4 for _ in range(4)]
    for N in (range(1, n_terms + 1) if orders is None else orders):
        a = xi ** N
        weight = -xi / (6 * mp.pi ** 2 * D) * (b * xi / (2 * D)) ** (N - 1)
        denominator = (b * (1 - mp.cos(theta)) * (1 + a) - T * (1 - a)) ** 3
        for m in rows:
            for n in range(4):
                out[m][n] += weight * (
                    (1 + 4 * a + a ** 2) * k[m] * k[n]
                    + (2 + 4 * a) * tan * (k[m] * rh[n] + rh[m] * k[n])
                    + 4 * tan ** 2 * rh[m] * rh[n]
                ) / denominator
    origin = _series_origin(orders, n_terms)
    src = KIM_THORNE_1991.ref(f"eq. 51 with eq. 52 ({variant})") + ("" if origin is Origin.DERIVED else ", series truncated at kt_series_tol")
    return [[Quantity(out[m][n], origin, src) for n in range(4)] for m in range(4)]


def klinkhammer_kk_on_axis(delta_t, X, D, b, xi, n_terms) -> Quantity:
    """Coefficient of k_a k_b from Klinkhammer (1992) eq. 3.8, specialized to this wormhole on the axis.

    Inputs from Kim & Thorne: Delta_n^(1/2) / zeta_n = D / [X (D - X)] (b / 2D)^n (their eq. 64 over
    zeta_n), and e^h = 1 / xi (their eq. 33). delta_t is the observer's proper time to H_n, X - T on the axis.
    """
    eh = 1 / xi
    total = mp.mpf(0)
    for n in range(1, n_terms + 1):
        amplitude = D / (X * (D - X)) * (b / (2 * D)) ** n
        total += amplitude * (1 + 4 * eh ** n + eh ** (2 * n)) / (eh ** n - 1) ** 3
    from .sources import KLINKHAMMER_1992

    origin = _series_origin(None, n_terms)
    return Quantity(-total / (12 * mp.pi ** 2 * delta_t ** 3), origin,
                    f"{KLINKHAMMER_1992.ref('eq. 3.8')} with {KIM_THORNE_1991.ref('eqs. 33, 64')}"
                    + ("" if origin is Origin.DERIVED else ", series truncated at kt_series_tol"))


# F2: geometric-optics amplitudes (eqs. 24, 47, 48, 57, 64)


def psi_throat(N, D, b):
    """psi_N = (1 / D) (b / 2D)^(N-1), eq. 47."""
    return (b / (2 * D)) ** (N - 1) / D


def vleck_throat(N, D, b, xi):
    """Delta_N^(1/2) = (zeta_N / D) (b / 2D)^(N-1), eq. 48, with zeta_N of eq. 36."""
    return zeta(N, D, xi) / D * (b / (2 * D)) ** (N - 1)


def zeta_between(N, X, D, xi):
    """zeta_N = D - X + sum_(n=1)^(N-1) D xi^n + X xi^N, eq. 57, closed form."""
    return (1 - xi ** N) / (1 - xi) * D - (1 - xi ** N) * X


def zeta_between_sum(N, X, D, xi):
    """The sum form of eq. 57."""
    return D - X + sum(D * xi ** n for n in range(1, N)) + X * xi ** N


def vleck_between(N, X, D, b, xi):
    """Delta_N^(1/2) = zeta_N D / [X (D - X)] (b / 2D)^N, eq. 64."""
    return zeta_between(N, X, D, xi) * D / (X * (D - X)) * (b / (2 * D)) ** N


# F5: polarized hypersurfaces and the Cauchy horizon (eqs. 43, 44a, 62, 63)


def time_HN_throat(N, theta, b, xi):
    """T_HN(theta) = [(1 + xi^N) / (1 - xi^N)] b (1 - cos theta), eq. 43."""
    return hypersurface_time_throat(N, theta, b, xi)


def time_Hminus_throat(theta, b):
    """T_H- = b (1 - cos theta), eq. 44a."""
    return b * (1 - mp.cos(theta))


def time_HN_between(N, X, rho, D, xi):
    """T_HN = X + rho^2 / [2X (D - X)] (D / (1 - xi^N) - X), eq. 62."""
    return X + rho ** 2 / (2 * X * (D - X)) * (D / (1 - xi ** N) - X)


def time_Hminus_between(X, rho):
    """T_H- = X + rho^2 / (2X), eq. 63."""
    return X + rho ** 2 / (2 * X)


# F6: four-momentum (eqs. 53, 54, 68, 69), natural units l_P = m_P = 1


def momentum_throat_N(N, delta_T, a, D, b, xi) -> Quantity:
    """k^mu coefficient of P^mu_N through the throat up to delta_T before H_N, eq. 53."""
    x = xi ** N
    value = -2 * b / (3 * mp.pi) / delta_T ** 2 * mp.log(b / a) * (b * xi / (2 * D)) ** N * (1 + 2 * x) / (1 - x) ** 3
    return Quantity(value, Origin.DERIVED, KIM_THORNE_1991.ref("eq. 53"))


def momentum_horizon(delta_T, D, b, xi, n_terms) -> Quantity:
    """k^mu coefficient of P^mu_- up to delta_T before H_-, eqs. 54 and 68 (the same expression)."""
    total = sum(
        (b / (2 * D)) ** N * (1 + 4 * xi ** N + xi ** (2 * N)) / (1 - xi ** N) ** 2 for N in range(1, n_terms + 1)
    )
    origin = _series_origin(None, n_terms)
    return Quantity(-total / (6 * mp.pi * delta_T), origin, KIM_THORNE_1991.ref("eqs. 54, 68")
                    + ("" if origin is Origin.DERIVED else f", sum truncated at N = {n_terms}"))


def momentum_per_area_N(N, X, D, b, xi) -> Quantity:
    """k^mu coefficient of (dP/dA)_N up to one Planck time before H_N, eq. 69."""
    x = xi ** N
    value = -D / (12 * mp.pi ** 2 * X * (D - X)) * (b * xi / (2 * D)) ** N * (1 + 4 * x + x ** 2) / (1 - x) ** 3
    return Quantity(value, Origin.DERIVED, KIM_THORNE_1991.ref("eq. 69"))
