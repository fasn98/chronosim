"""Independent numerics for the wormhole spacetime of Kim & Thorne (1991).

The covering space is built here directly from their eqs. 25 and 28, as repeated
Poincare maps of the T-X plane, without their eqs. 29-38. Everything here is
COMPUTED or a building block of a COMPUTED result.
"""

from __future__ import annotations

import mpmath as mp

from .config import DPS
from .kim_thorne import hadamard_flat
from .pointsplit import DELTA, STEP
from .provenance import Origin, Quantity
from .sources import KIM_THORNE_1991

METRIC = (-1, 1, 1, 1)


def _shift(P, i, s):
    P = list(P)
    P[i] = P[i] + s
    return P


def _eq4_at(F, x, xp, h):
    """Eq. 4b applied to F(x, x') at a split pair, by central differences."""
    mixed = [[(F(_shift(x, a, h), _shift(xp, c, h)) - F(_shift(x, a, h), _shift(xp, c, -h))
               - F(_shift(x, a, -h), _shift(xp, c, h)) + F(_shift(x, a, -h), _shift(xp, c, -h))) / (4 * h * h)
              for c in range(4)] for a in range(4)]
    same = [[(F(_shift(_shift(x, a, h), c, h), xp) - F(_shift(_shift(x, a, h), c, -h), xp)
              - F(_shift(_shift(x, a, -h), c, h), xp) + F(_shift(_shift(x, a, -h), c, -h), xp)) / (4 * h * h)
             for c in range(4)] for a in range(4)]
    contraction = sum(METRIC[c] * mixed[c][c] for c in range(4))
    return [[mp.mpf(2) / 3 * mixed[a][c] - mp.mpf(1) / 3 * same[a][c]
             - (mp.mpf(1) / 6 * METRIC[a] * contraction if a == c else 0)
             for c in range(4)] for a in range(4)]


def eq4_stress(F, X, delta=None, h=None, dps=None):
    """Cartesian T_ab of Kim & Thorne eqs. 4a and 4b, exactly as printed, for a regularized G.

    Written from those equations and not from ``pointsplit.stress_tensor``, so that
    the ratio between the two conventions is measured here. The split is along y,
    where ``stress_tensor`` splits along z, and +delta and -delta are averaged.
    """
    delta = DELTA.value if delta is None else mp.mpf(delta)
    h = STEP.value if h is None else mp.mpf(h)
    with mp.workdps(DPS.value if dps is None else dps):
        X = [mp.mpf(v) for v in X]
        plus = _eq4_at(F, X, _shift(X, 2, delta), h)
        minus = _eq4_at(F, X, _shift(X, 2, -delta), h)
        return [[(plus[a][c] + minus[a][c]) / 2 for c in range(4)] for a in range(4)]


def half_square_separation(P, Q):
    """sigma_0 of Kim & Thorne eq. 3: half the square of the separation (MTW signs)."""
    dt, dx, dy, dz = (mp.mpf(a) - mp.mpf(c) for a, c in zip(P, Q))
    return (-dt ** 2 + dx ** 2 + dy ** 2 + dz ** 2) / 2


def regularized_with_eq3(G):
    """G_reg = G - 1 / (4 pi^2 sigma_0), Kim & Thorne eq. 3."""
    return lambda P, Q: G(P, Q) - hadamard_flat(half_square_separation(P, Q))


def throat_map(event, D, beta):
    """Carry copy-n coordinates to copy n + 1: the Poincare map sending the first world
    line to the second at equal proper time (eqs. 25a, 25b and 28).

    (T = tau, X = 0) goes to (T = D + gamma tau, X = D - beta gamma tau).
    """
    T, X = event
    gamma = 1 / mp.sqrt(1 - beta ** 2)
    return D + gamma * (T - beta * X), D + gamma * (X - beta * T)


def copy_event(event, n, D, beta):
    """Copy n of an event, in the copy 0 coordinates, by iterating ``throat_map``."""
    for _ in range(n):
        event = throat_map(event, D, beta)
    return event


def exact_interval(N, x, xp, D, beta):
    """Half the square of the straight line from x' in copy 0 to copy N of x, with MTW signs."""
    T, X = copy_event(x, N, D, beta)
    return (-(T - xp[0]) ** 2 + (X - xp[1]) ** 2) / 2


def exact_sigma_N(N, x, xp, D, beta) -> Quantity:
    """Geodetic interval on the symmetry axis, without the approximation of eq. 35."""
    with mp.workdps(DPS.value):
        value = exact_interval(N, x, xp, mp.mpf(D), mp.mpf(beta))
    return Quantity(value, Origin.COMPUTED, f"covering space built from {KIM_THORNE_1991.ref('eqs. 25, 28')}")
