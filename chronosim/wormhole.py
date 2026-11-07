"""Independent numerics for the wormhole spacetime of Kim & Thorne (1991).

The covering space is built here directly from their eqs. 25 and 28, as repeated
Poincare maps of the T-X plane, without their eqs. 29-38. Everything here is
COMPUTED or a building block of a COMPUTED result.
"""

from __future__ import annotations

import mpmath as mp

from .config import DPS
from . import kim_thorne as kt
from .kim_thorne import doppler_xi, hadamard_between, hadamard_flat, hadamard_throat, series_terms
from .pointsplit import DELTA, STEP, stress_tensor
from .provenance import Origin, Quantity, heuristic
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


SPLIT_RATIO = heuristic("kt_split_ratio", mp.mpf("1e-10"), "point-splitting distance / distance to the horizon, Kim & Thorne runs")
STEP_RATIO = heuristic("kt_step_ratio", mp.mpf("1e-11"), "finite-difference step / distance to the horizon, Kim & Thorne runs")


def stress_between_numeric(T, X, rho, D, b, beta, length, tol=None, split_ratio=None, step_ratio=None) -> dict:
    """T^mu^nu between the mouths from eq. 66, without using eq. 67, with two operators.

    ``eq4a`` is ``eq4_stress`` (eq. 4a as printed, split along y); ``with_half`` is
    ``pointsplit.stress_tensor`` (the same operator with the factor 1/2, split along z).
    ``length`` is the smallest length scale, the distance to the horizon, and sets
    the split and the step. Contravariant components in (T, X) are returned as
    TT, TX, XX, plus the largest remaining component.
    """
    split = (SPLIT_RATIO.value if split_ratio is None else mp.mpf(split_ratio)) * length
    step = (STEP_RATIO.value if step_ratio is None else mp.mpf(step_ratio)) * length
    with mp.workdps(DPS.value):
        xi = doppler_xi(mp.mpf(beta))
        n_terms = series_terms(D, b, xi, tol)
        G = lambda P, Q: hadamard_between(P, Q, D, b, xi, n_terms)
        point = [T, X, rho, mp.mpf(0)]
        tensors = {
            "eq4a": eq4_stress(G, point, delta=split, h=step),
            "with_half": stress_tensor(G, point, delta=split, h=step),
        }
        out = {"n_terms": n_terms}
        for name, L in tensors.items():
            up = [[L[a][c] * METRIC[a] * METRIC[c] for c in range(4)] for a in range(4)]
            rest = max(abs(up[a][c]) for a in range(4) for c in range(4) if a > 1 or c > 1)
            out[name] = {
                "TT": Quantity(up[0][0], Origin.COMPUTED, f"{KIM_THORNE_1991.ref('eq. 66')} + {name} operator"),
                "TX": Quantity(up[0][1], Origin.COMPUTED, f"{KIM_THORNE_1991.ref('eq. 66')} + {name} operator"),
                "XX": Quantity(up[1][1], Origin.COMPUTED, f"{KIM_THORNE_1991.ref('eq. 66')} + {name} operator"),
                "rest": Quantity(rest, Origin.COMPUTED, f"{KIM_THORNE_1991.ref('eq. 66')} + {name} operator"),
            }
    return out


def stress_throat_numeric(point, D, b, beta, length, tol=None, split_ratio=None, step_ratio=None) -> dict:
    """T^mu^nu near the throat from eq. 50, without using eq. 51, with both operators."""
    split = (SPLIT_RATIO.value if split_ratio is None else mp.mpf(split_ratio)) * length
    step = (STEP_RATIO.value if step_ratio is None else mp.mpf(step_ratio)) * length
    with mp.workdps(DPS.value):
        xi = doppler_xi(mp.mpf(beta))
        n_terms = series_terms(D, b, xi, tol)
        G = lambda P, Q: hadamard_throat(P, Q, D, b, xi, n_terms)
        point = [mp.mpf(v) for v in point]
        out = {"n_terms": n_terms}
        for name, L in (("eq4a", eq4_stress(G, point, delta=split, h=step)),
                        ("with_half", stress_tensor(G, point, delta=split, h=step))):
            src = f"{KIM_THORNE_1991.ref('eq. 50')} + {name} operator"
            out[name] = [[Quantity(L[a][c] * METRIC[a] * METRIC[c], Origin.COMPUTED, src) for c in range(4)]
                         for a in range(4)]
    return out


# F2: ray-transfer (ABCD) matrices for the defocusing, independent of eqs. 47, 48 and 64.
# The wormhole is the thin diverging lens of focal length b / 2 that the text describes
# (p. 3942). A point source has amplitude 1 / B, with B the upper-right entry.


def _matmul(A, B):
    return [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
            [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]]


def _propagate(L):
    return [[mp.mpf(1), L], [mp.mpf(0), mp.mpf(1)]]


def _diverging_lens(f):
    return [[mp.mpf(1), mp.mpf(0)], [1 / f, mp.mpf(1)]]


def ray_amplitude_throat(N, D, b) -> Quantity:
    """Amplitude after a distance D and N - 1 more (lens, distance D) legs, for a source at mouth 1."""
    with mp.workdps(DPS.value):
        M = _propagate(D)
        for _ in range(N - 1):
            M = _matmul(_propagate(D), _matmul(_diverging_lens(b / 2), M))
    return Quantity(1 / M[0][1], Origin.COMPUTED, "ABCD matrices, lens of focal length b/2 (Kim & Thorne p. 3942)")


def ray_amplitude_between(N, X, D, b) -> Quantity:
    """Amplitude for a source between the mouths: D - X, then N (lens, D) legs with the last of length X."""
    with mp.workdps(DPS.value):
        M = _propagate(D - X)
        for n in range(N):
            M = _matmul(_propagate(X if n == N - 1 else D), _matmul(_diverging_lens(b / 2), M))
    return Quantity(1 / M[0][1], Origin.COMPUTED, "ABCD matrices, lens of focal length b/2 (Kim & Thorne p. 3942)")


def exact_lapse_on_horizon(N, X, D, beta):
    """V_N(x) - V_0(x) for x on the horizon (T = X), from the covering space of eqs. 25 and 28."""
    with mp.workdps(DPS.value):
        xi = doppler_xi(mp.mpf(beta))
        x = (mp.mpf(X), mp.mpf(X))
        _, V = kt.null_coords(kt.to_covering(copy_event(x, N, D, mp.mpf(beta)), D, xi))
        _, V0 = kt.null_coords(kt.to_covering(x, D, xi))
        return V - V0


# F5: hypersurfaces as roots of the collapsed lambda_N, independent of eqs. 43 and 62.


def root_HN_throat(N, theta, b, beta) -> Quantity:
    """T where lambda_N(x, x) of eq. 46 vanishes on the throat point of eq. 42."""
    with mp.workdps(DPS.value):
        xi = doppler_xi(mp.mpf(beta))
        f = lambda T: kt.lambda_throat(N, kt.throat_point(T, theta, 0, b), kt.throat_point(T, theta, 0, b), b, xi)
        value = mp.findroot(f, b)
    return Quantity(value, Origin.COMPUTED, f"root of {KIM_THORNE_1991.ref('eq. 46')}")


def root_HN_between(N, X, rho, D, beta) -> Quantity:
    """T where lambda_N(x, x) of eq. 61 vanishes."""
    with mp.workdps(DPS.value):
        xi = doppler_xi(mp.mpf(beta))
        f = lambda T: kt.lambda_between(N, [T, X, rho, 0], [T, X, rho, 0], D, xi)
        value = mp.findroot(f, X)
    return Quantity(value, Origin.COMPUTED, f"root of {KIM_THORNE_1991.ref('eq. 61')}")


# F6: four-momentum by numerical integration of eqs. 51 and 67, independent of eqs. 53, 54, 68, 69.

QUAD_DPS = heuristic("kt_quad_dps", 15, "working precision of the F6 integrals")
THROAT_THICKNESS = heuristic("kt_throat_thickness", mp.mpf("1e-6"), "a / b, angular cutoff tan(theta) < b / a for eqs. 53-54")


def _throat_flux(T, theta, D, b, xi, n_terms, orders):
    """T^T^nu n_nu on the throat, from eq. 51 with the transverse rho_hat, at phi = 0."""
    E = kt.stress_throat(T, theta, mp.mpf(0), D, b, xi, n_terms, "transverse", orders, rows=(0,))
    n = (mp.mpf(0), mp.cos(theta), mp.sin(theta), mp.mpf(0))
    return sum(E[0][c].value * n[c] for c in range(4))


def momentum_throat_numeric(delta_T, D, b, beta, upper, orders=None, thickness=None) -> Quantity:
    """k^T component of the integral of eq. 51 over the inner face and T, up to ``upper(theta) - delta_T``.

    The integrand does not depend on phi, so the phi integral is 2 pi.
    """
    a_over_b = THROAT_THICKNESS.value if thickness is None else mp.mpf(thickness)
    with mp.workdps(QUAD_DPS.value):
        xi = doppler_xi(mp.mpf(beta))
        n_terms = series_terms(D, b, xi)
        theta_max = mp.atan(1 / a_over_b)

        def over_T(theta):
            top = upper(theta) - delta_T
            return mp.quad(lambda s: _throat_flux(top - s, theta, D, b, xi, n_terms, orders), [0, delta_T, 100 * delta_T, mp.inf])

        breaks = [0, mp.sqrt(delta_T / b) / 10, mp.sqrt(delta_T / b), mp.mpf("0.5"), mp.pi / 2 - mp.mpf("0.1"), theta_max]
        value = 2 * mp.pi * b ** 2 * mp.quad(lambda th: over_T(th) * mp.sin(th), breaks)
    return Quantity(value, Origin.COMPUTED, f"quadrature of {KIM_THORNE_1991.ref('eq. 51')}")


def momentum_plane_numeric(delta_T, X, D, b, beta, orders=None) -> Quantity:
    """k^T component of the integral of T^T^X from eq. 67 over the plane and T, up to T_H- - delta_T (eq. 63)."""
    with mp.workdps(QUAD_DPS.value):
        xi = doppler_xi(mp.mpf(beta))
        n_terms = series_terms(D, b, xi)

        def over_T(rho):
            top = kt.time_Hminus_between(X, rho) - delta_T
            return mp.quad(lambda s: kt.stress_between_kk(top - s, X, rho, D, b, xi, n_terms, orders).value,
                           [0, delta_T, 100 * delta_T, mp.inf])

        scale = mp.sqrt(delta_T * X)
        value = mp.quad(lambda r: over_T(r) * 2 * mp.pi * r, [0, scale, 10 * scale, 100 * scale, mp.inf])
    return Quantity(value, Origin.COMPUTED, f"quadrature of {KIM_THORNE_1991.ref('eq. 67')}")


def momentum_per_area_numeric(N, X, rho, D, b, beta, cutoff=1) -> Quantity:
    """k^T component of the integral of the order-N term of eq. 67 in T, up to T_HN - cutoff (eq. 62)."""
    with mp.workdps(QUAD_DPS.value):
        xi = doppler_xi(mp.mpf(beta))
        top = kt.time_HN_between(N, X, rho, D, xi) - cutoff
        value = mp.quad(lambda s: kt.stress_between_kk(top - s, X, rho, D, b, xi, N, [N]).value,
                        [0, cutoff, 100 * cutoff, mp.inf])
    return Quantity(value, Origin.COMPUTED, f"quadrature of {KIM_THORNE_1991.ref('eq. 67')}")
