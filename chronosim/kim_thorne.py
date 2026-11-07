"""Formulas transcribed from Kim & Thorne (1991). See docs/equation_transcription.md.

Conventions (their Sec. I): c = G = hbar = 1, signs of Misner, Thorne & Wheeler.
Their xi is the inverse Doppler blueshift of eq. 30, not the Rindler coordinate
of ``misner.py``. Events in the T-X plane are (T, X).
"""

from __future__ import annotations

import mpmath as mp

from .provenance import Origin, Quantity
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
