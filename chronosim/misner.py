"""Misner space: Rindler coordinates, the image-sum Hadamard function, and the
closed-form renormalized stress-energy of Li & Gott (1998).

Conventions (Li & Gott 1998, eqs. 1-3): signature (-,+,+,+); Cartesian
X = (t, x, y, z); in the right Rindler wedge R (x > |t|)
    t = xi sinh(eta),  x = xi cosh(eta),  ds^2 = -xi^2 deta^2 + dxi^2 + dy^2 + dz^2.
Misner space identifies eta ~ eta + n b (a boost by rapidity b). In R the
Killing vector d/deta is timelike, so these identifications create closed
timelike curves; the chronology horizon is xi -> 0 (the light cone x = |t|).
"""

from __future__ import annotations

import mpmath as mp

from .config import DPS
from .provenance import Origin, Quantity
from .sources import LI_GOTT_1998


def to_rindler(t, x):
    """(t, x) -> (xi, eta) in region R. Li & Gott eq. 2 inverted."""
    t, x = mp.mpf(t), mp.mpf(x)
    if not x > abs(t):
        raise ValueError("point is outside the right Rindler wedge (need x > |t|)")
    return mp.sqrt((x - t) * (x + t)), mp.atanh(t / x)


def _gamma(xi, xip, dy, dz):
    """Return (gamma^2, gamma / sinh gamma) with cosh(gamma) from eq. 7.

    cosh(gamma) - 1 = [(xi - xi')^2 + dy^2 + dz^2] / (2 xi xi') is computed
    directly, and gamma = 2 asinh(sqrt(u/2)), which is exact and avoids the
    cancellation of acosh near 1.
    """
    u = ((xi - xip) ** 2 + dy ** 2 + dz ** 2) / (2 * xi * xip)
    if u == 0:
        return mp.mpf(0), mp.mpf(1)
    g = 2 * mp.asinh(mp.sqrt(u / 2))
    return g * g, g / mp.sinh(g)


def hadamard_image(n: int, X, Xp, b):
    """n-th term of the Misner-space Hadamard function, Li & Gott eq. 8.

    The n = 0 term is the Rindler-vacuum Hadamard function (their eq. 6).
    """
    t, x, y, z = X
    tp, xp, yp, zp = Xp
    xi, eta = to_rindler(t, x)
    xip, etap = to_rindler(tp, xp)
    g2, ratio = _gamma(xi, xip, y - yp, z - zp)
    return ratio / (2 * mp.pi ** 2 * xi * xip * (-(eta - etap + n * b) ** 2 + g2))


def hadamard_minkowski(X, Xp):
    """Minkowski-vacuum Hadamard function, Li & Gott eq. 9."""
    dt, dx, dy, dz = (mp.mpf(a) - mp.mpf(c) for a, c in zip(X, Xp))
    return 1 / (2 * mp.pi ** 2 * (-dt ** 2 + dx ** 2 + dy ** 2 + dz ** 2))


def boost(rapidity, X):
    """Boost in the (t, x) plane, eta -> eta + rapidity. This is the Misner identification."""
    t, x, y, z = X
    c, s = mp.cosh(rapidity), mp.sinh(rapidity)
    return [t * c + x * s, t * s + x * c, y, z]


def hadamard_minkowski_image(n: int, X, Xp, b):
    """n-th image of the Minkowski vacuum, G_M(X, Lambda^n X'), with Lambda the boost of rapidity b.

    Summing this over n gives the Minkowski vacuum with images. The n = 0 term is
    what the regularization subtracts. G_M is Li & Gott eq. 9.
    """
    return hadamard_minkowski(X, boost(n * b, Xp))


def stress_prefactor(xi, b):
    """C = [(2 pi / b)^4 - 1] / (1440 pi^2 xi^4), the scalar factor of eq. 12."""
    with mp.workdps(DPS.value):
        xi, b = mp.mpf(xi), mp.mpf(b)
        return ((2 * mp.pi / b) ** 4 - 1) / (1440 * mp.pi ** 2 * xi ** 4)


def closed_form_stress(xi, b) -> dict:
    """Li & Gott eq. 12 in the orthonormal Rindler frame (e_eta, e_xi, e_y, e_z).

    <T_ab> = C (4 xi^2 deta_a deta_b + g_ab). With deta_a deta^a = -1/xi^2 the
    frame components are: energy density rho = 3C, pressures p_xi = p_y = p_z = C.
    """
    C = stress_prefactor(xi, b)
    src = LI_GOTT_1998.ref("eq. 12")
    with mp.workdps(DPS.value):
        rho = 3 * C
    return {
        "rho": Quantity(rho, Origin.DERIVED, src),
        "p_xi": Quantity(C, Origin.DERIVED, src),
        "p_y": Quantity(C, Origin.DERIVED, src),
        "p_z": Quantity(C, Origin.DERIVED, src),
    }


def self_consistent_b() -> Quantity:
    """b = 2 pi, the value for which Li & Gott find <T_ab> = 0 everywhere."""
    with mp.workdps(DPS.value):
        return Quantity(+(2 * mp.pi), Origin.CITED, LI_GOTT_1998.ref("eq. 13: b = 2 pi gives <T_ab> = 0"))
