"""Renormalized stress-energy by point splitting, Li & Gott (1998) eq. 11:

    <T_ab> = 1/2 lim_{X'->X} [ 2/3 d_a d_b'  -  1/3 d_a d_b  -  1/6 g_ab d_c d^c' ] G_reg(X, X')

In flat space with Cartesian coordinates covariant derivatives are partial
derivatives and the parallel propagator is the identity. Derivatives are taken
by central finite differences in arbitrary precision; the limit X' -> X is
approached by splitting along z by +delta and -delta and averaging, which
cancels the odd terms. Everything here is a numerical approximation, so every
output is COMPUTED.
"""

from __future__ import annotations

import mpmath as mp

from .misner import hadamard_image, hadamard_minkowski, hadamard_minkowski_image, to_rindler
from .config import DPS
from .provenance import Origin, Quantity, heuristic
from .sources import LI_GOTT_1998

METRIC = (-1, 1, 1, 1)

DELTA = heuristic("delta", mp.mpf("1e-10"), "point-splitting distance along z, units of length")
STEP = heuristic("h", mp.mpf("1e-11"), "finite-difference step (must be < delta), units of length")


def _shift(P, i, s):
    P = list(P)
    P[i] = P[i] + s
    return P


def _tensor_at(F, X, Xp, h):
    """Apply the eq. 11 operator to F at the split pair (X, X'). Returns 4x4 list."""
    f0 = F(X, Xp)
    A = [[None] * 4 for _ in range(4)]  # d_a d_b'
    B = [[None] * 4 for _ in range(4)]  # d_a d_b (first argument)
    for a in range(4):
        for c in range(4):
            pp = F(_shift(X, a, h), _shift(Xp, c, h))
            pm = F(_shift(X, a, h), _shift(Xp, c, -h))
            mp_ = F(_shift(X, a, -h), _shift(Xp, c, h))
            mm = F(_shift(X, a, -h), _shift(Xp, c, -h))
            A[a][c] = (pp - pm - mp_ + mm) / (4 * h * h)
        for c in range(a, 4):
            if a == c:
                v = (F(_shift(X, a, h), Xp) - 2 * f0 + F(_shift(X, a, -h), Xp)) / (h * h)
            else:
                XA = _shift(X, a, h)
                XB = _shift(X, a, -h)
                v = (
                    F(_shift(XA, c, h), Xp) - F(_shift(XA, c, -h), Xp)
                    - F(_shift(XB, c, h), Xp) + F(_shift(XB, c, -h), Xp)
                ) / (4 * h * h)
            B[a][c] = B[c][a] = v
    trA = sum(METRIC[c] * A[c][c] for c in range(4))  # g^{cd} d_c d_d'
    T = [[None] * 4 for _ in range(4)]
    for a in range(4):
        for c in range(4):
            g = METRIC[a] if a == c else 0
            T[a][c] = (mp.mpf(2) / 3 * A[a][c] - mp.mpf(1) / 3 * B[a][c] - mp.mpf(1) / 6 * g * trA) / 2
    return T


def stress_tensor(F, X, delta=None, h=None, dps=None):
    """Cartesian <T_ab> at X for a regularized Hadamard function F(X, X')."""
    delta = DELTA.value if delta is None else mp.mpf(delta)
    h = STEP.value if h is None else mp.mpf(h)
    if not h < delta:
        raise ValueError("finite-difference step must be smaller than the split distance")
    with mp.workdps(DPS.value if dps is None else dps):
        X = [mp.mpf(v) for v in X]
        Tp = _tensor_at(F, X, _shift(X, 3, delta), h)
        Tm = _tensor_at(F, X, _shift(X, 3, -delta), h)
        return [[(Tp[a][c] + Tm[a][c]) / 2 for c in range(4)] for a in range(4)]


def rindler_part(X, **kw):
    """n = 0 image minus the Minkowski vacuum (the Rindler-vacuum contribution)."""
    F = lambda P, Q: hadamard_image(0, P, Q, 1) - hadamard_minkowski(P, Q)  # b irrelevant for n=0
    return stress_tensor(F, X, **kw)


N_TERMS = heuristic("n_terms", 24, "image pairs summed explicitly before the Kummer tail")


def image_pair(X, b, n, hadamard=hadamard_image, **kw):
    """T(n) + T(-n) for the images of ``hadamard`` (a 4x4 list)."""
    Tn = stress_tensor(lambda P, Q: hadamard(n, P, Q, b), X, **kw)
    Tmn = stress_tensor(lambda P, Q: hadamard(-n, P, Q, b), X, **kw)
    return [[Tn[a][c] + Tmn[a][c] for c in range(4)] for a in range(4)]


def image_pairs(X, b, n_max, **kw):
    """T(n) + T(-n) for n = 1..n_max (each a 4x4 list)."""
    b = mp.mpf(b)
    return [image_pair(X, b, n, **kw) for n in range(1, n_max + 1)]


def image_part(X, b, n_terms=None, tail=True, diagnostics=None, **kw):
    """Sum over n != 0 of the eq. 8 images, sum_{n>=1} [T(n) + T(-n)].

    The first ``n_terms`` pairs are summed explicitly. The remainder uses
    Kummer's acceleration: assuming pair(n) ~ K / n^4 for large n, the tail is
    K_N * (zeta(4) - sum_{n<=N} n^-4) with K_N = N^4 pair(N). The assumption
    is checked, not trusted: ``diagnostics['decay_spread']`` reports the
    largest relative spread of n^4 pair(n) over the last half of the explicit
    terms (it should be at the level of the finite-difference error).
    """
    N = N_TERMS.value if n_terms is None else n_terms
    with mp.workdps(DPS.value if kw.get("dps") is None else kw["dps"]):
        pairs = image_pairs(X, b, N, **kw)
        S = [[sum(p[a][c] for p in pairs) for c in range(4)] for a in range(4)]
        if not tail:
            return S
        rest = mp.zeta(4) - sum(mp.mpf(n) ** -4 for n in range(1, N + 1))
        K = [[pairs[-1][a][c] * N ** 4 for c in range(4)] for a in range(4)]
        if diagnostics is not None:
            spread = mp.mpf(0)
            for a in range(4):
                ref = K[a][a]
                if ref == 0:
                    continue
                for n in range(N // 2, N + 1):
                    spread = max(spread, abs(pairs[n - 1][a][a] * n ** 4 - ref) / abs(ref))
            diagnostics["decay_spread"] = spread
            diagnostics["tail_fraction"] = max(
                abs(K[a][a] * rest / (S[a][a] + K[a][a] * rest)) for a in range(4) if S[a][a] != 0
            )
        return [[S[a][c] + K[a][c] * rest for c in range(4)] for a in range(4)]


MINK_TOL = heuristic("mink_tol", mp.mpf("1e-25"), "Minkowski image sum stops when |pair_00| <= mink_tol |sum_00|")


def minkowski_image_part(X, b, tol=None, diagnostics=None, **kw):
    """Sum over n != 0 of the Minkowski-vacuum images, sum_{n>=1} [T(n) + T(-n)].

    The decay law was measured, not assumed, and it is not the n^-4 of
    ``image_part``: pair(n) looks like n^-4 only while n b << 1, then the ratio
    r_n = pair(n) / pair(n-1) settles to a constant with -ln(r_n) / b -> 1. So the
    sum runs until the last pair is below ``tol`` of the total and the remainder
    is the geometric tail pair(N) r / (1 - r), with r measured from the last two
    pairs. A fixed N would not do: b sets how many terms matter, and by n b ~ 140
    at dps = 80 cancellation in sigma^2 has turned the terms into noise.

    ``diagnostics`` receives ``decay_rate_over_b`` (-ln(r_N) / b), ``power_law_spread``
    (relative spread of n^4 pair(n) over the last half of the terms, large when
    n^-4 fails), ``tail_fraction`` and ``n_used``.
    """
    tol = MINK_TOL.value if tol is None else mp.mpf(tol)
    with mp.workdps(DPS.value if kw.get("dps") is None else kw["dps"]):
        b = mp.mpf(b)
        if not b > 0:
            raise ValueError("the boost rapidity b must be positive")
        pairs = [image_pair(X, b, 1, hadamard=hadamard_minkowski_image, **kw)]
        S = pairs[0]
        while True:
            pairs.append(image_pair(X, b, len(pairs) + 1, hadamard=hadamard_minkowski_image, **kw))
            S = [[S[a][c] + pairs[-1][a][c] for c in range(4)] for a in range(4)]
            if abs(pairs[-1][0][0]) <= tol * abs(S[0][0]):
                break
        N = len(pairs)
        r = pairs[-1][0][0] / pairs[-2][0][0]
        tail = [[pairs[-1][a][c] * r / (1 - r) for c in range(4)] for a in range(4)]
        if diagnostics is not None:
            weighted = [pairs[n - 1][0][0] * n ** 4 for n in range(N // 2 + 1, N + 1)]
            diagnostics["decay_rate_over_b"] = -mp.log(abs(r)) / b
            diagnostics["power_law_spread"] = (max(weighted) - min(weighted)) / abs(weighted[-1])
            diagnostics["tail_fraction"] = abs(tail[0][0] / (S[0][0] + tail[0][0]))
            diagnostics["n_used"] = mp.mpf(N)
        return [[S[a][c] + tail[a][c] for c in range(4)] for a in range(4)]


def minkowski_images_stress(xi, b, **kw) -> dict:
    """Numerical <T_ab> of the Minkowski vacuum with images, on the t = 0 slice at
    Rindler radius xi (see ``misner_stress``). All outputs are COMPUTED.
    """
    with mp.workdps(DPS.value if kw.get("dps") is None else kw["dps"]):
        X = [0, mp.mpf(xi), 0, 0]
        to_rindler(X[0], X[1])  # validates the point
        diag = {}
        T = minkowski_image_part(X, b, diagnostics=diag, **kw)
        trace = -T[0][0] + T[1][1] + T[2][2] + T[3][3]
    src = f"Minkowski image sum (G_M: {LI_GOTT_1998.ref('eq. 9')}) + point splitting ({LI_GOTT_1998.ref('eq. 11')})"
    names = {"rho": 0, "p_xi": 1, "p_y": 2, "p_z": 3}
    out = {k: Quantity(T[i][i], Origin.COMPUTED, src) for k, i in names.items()}
    out["max_offdiag"] = Quantity(
        max(abs(T[a][c]) for a in range(4) for c in range(4) if a != c), Origin.COMPUTED, src
    )
    out["trace"] = Quantity(trace, Origin.COMPUTED, src)
    for k, v in diag.items():
        out[k] = Quantity(v, Origin.COMPUTED, "Minkowski image-sum convergence diagnostic")
    return out


def misner_stress(xi, b, **kw) -> dict:
    """Numerical <T_ab> on the t = 0 slice at Rindler radius xi (so eta = 0 and the
    Cartesian frame coincides with the orthonormal Rindler frame).

    Returns frame components as COMPUTED quantities, plus the largest
    off-diagonal component as a symmetry check.
    """
    with mp.workdps(DPS.value if kw.get("dps") is None else kw["dps"]):
        X = [0, mp.mpf(xi), 0, 0]
        to_rindler(X[0], X[1])  # validates the point
        diag = {}
        R = rindler_part(X, **kw)
        I = image_part(X, b, diagnostics=diag, **kw)
        T = [[R[a][c] + I[a][c] for c in range(4)] for a in range(4)]
        trace = -T[0][0] + T[1][1] + T[2][2] + T[3][3]
    src = f"image sum ({LI_GOTT_1998.ref('eq. 8')}) + point splitting ({LI_GOTT_1998.ref('eq. 11')})"
    names = {"rho": 0, "p_xi": 1, "p_y": 2, "p_z": 3}
    out = {k: Quantity(T[i][i], Origin.COMPUTED, src) for k, i in names.items()}
    out["max_offdiag"] = Quantity(
        max(abs(T[a][c]) for a in range(4) for c in range(4) if a != c), Origin.COMPUTED, src
    )
    out["trace"] = Quantity(trace, Origin.COMPUTED, src)
    for k, v in diag.items():
        out[k] = Quantity(v, Origin.COMPUTED, "image-sum convergence diagnostic")
    return out
