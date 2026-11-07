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

from .misner import hadamard_image, hadamard_minkowski, to_rindler
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


def image_pairs(X, b, n_max, **kw):
    """T(n) + T(-n) for n = 1..n_max (each a 4x4 list)."""
    b = mp.mpf(b)
    out = []
    for n in range(1, n_max + 1):
        Tn = stress_tensor(lambda P, Q, n=n: hadamard_image(n, P, Q, b), X, **kw)
        Tmn = stress_tensor(lambda P, Q, n=n: hadamard_image(-n, P, Q, b), X, **kw)
        out.append([[Tn[a][c] + Tmn[a][c] for c in range(4)] for a in range(4)])
    return out


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
