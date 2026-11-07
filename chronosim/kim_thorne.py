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


def stress_between_kk(T, X, rho, D, b, xi, n_terms) -> Quantity:
    """Coefficient of k^mu k^nu in T^mu^nu between the mouths, eq. 67, as printed."""
    total = mp.mpf(0)
    for N in range(1, n_terms + 1):
        a = xi ** N
        denominator = (X - T) * (1 - a) + rho ** 2 / (2 * X) + rho ** 2 / (2 * (D - X)) * a
        total += (b * xi / (2 * D)) ** N * (1 + 4 * a + a ** 2) / denominator ** 3
    return Quantity(-D / (6 * mp.pi ** 2 * X * (D - X)) * total, Origin.DERIVED, KIM_THORNE_1991.ref("eq. 67"))
