"""Arbiter that depends on no closed article: a thermal state in Minkowski space.

The Hadamard function of the state at inverse temperature beta is the image sum in
imaginary time derived from the KMS condition in docs/thermal_images.md. The n = 0
term is the Minkowski vacuum, which the regularization removes.
"""

from __future__ import annotations

import mpmath as mp

from .config import DPS
from .pointsplit import stress_tensor
from .provenance import Origin, Quantity, heuristic
from .wormhole import eq4_stress

N_TERMS = heuristic("thermal_n_terms", 24, "imaginary-time image pairs summed explicitly before the Kummer tail")
THERMAL_SPLIT_RATIO = heuristic("thermal_split_ratio", mp.mpf("1e-10"), "point-splitting distance / beta in the thermal arbiter")
THERMAL_STEP_RATIO = heuristic("thermal_step_ratio", mp.mpf("1e-11"), "finite-difference step / beta in the thermal arbiter")


MODE_CUTOFF = heuristic("thermal_mode_cutoff", 60, "mode integral stops at k = cutoff / beta, where n_k ~ e^-60")


def hadamard_image_pair(n, P, Q, beta):
    """G_M at t - t' + i n beta plus the same at -i n beta, a real number. n >= 1."""
    dt, dx, dy, dz = (mp.mpf(a) - mp.mpf(c) for a, c in zip(P, Q))
    r2 = dx ** 2 + dy ** 2 + dz ** 2
    image = 1 / (2 * mp.pi ** 2 * (r2 - mp.mpc(dt, n * beta) ** 2))
    return 2 * image.real


def hadamard_regularized(P, Q, beta, n_max):
    """Sum over the image pairs up to n_max, without tail."""
    return sum(hadamard_image_pair(n, P, Q, beta) for n in range(1, n_max + 1))


def hadamard_mode_integral(t, r, beta) -> Quantity:
    """G_reg from the thermal occupation numbers, (1 / 2 pi^2) int dk 2 n_k cos(k t) sin(k r) / r."""
    with mp.workdps(30):
        value = mp.quad(
            lambda k: 2 / mp.expm1(beta * k) * mp.cos(k * t) * mp.sin(k * r) / r,
            mp.linspace(0, MODE_CUTOFF.value / beta, 31),
        ) / (2 * mp.pi ** 2)
    return Quantity(value, Origin.COMPUTED, "mode integral with Bose occupation numbers, docs/thermal_images.md")


def _sum_operator(operator, X, beta, n_terms, split, step, diagnostics):
    pairs = [
        operator(lambda P, Q, n=n: hadamard_image_pair(n, P, Q, beta), X, delta=split, h=step)
        for n in range(1, n_terms + 1)
    ]
    S = [[sum(p[a][c] for p in pairs) for c in range(4)] for a in range(4)]
    rest = mp.zeta(4) - sum(mp.mpf(n) ** -4 for n in range(1, n_terms + 1))
    K = [[pairs[-1][a][c] * n_terms ** 4 for c in range(4)] for a in range(4)]
    if diagnostics is not None:
        weighted = [pairs[n - 1][0][0] * n ** 4 for n in range(n_terms // 2, n_terms + 1)]
        diagnostics["decay_spread"] = (max(weighted) - min(weighted)) / abs(weighted[-1])
        diagnostics["tail_fraction"] = abs(K[0][0] * rest / (S[0][0] + K[0][0] * rest))
    return [[S[a][c] + K[a][c] * rest for c in range(4)] for a in range(4)]


def thermal_stress(T, X=None, n_terms=None, split_ratio=None, step_ratio=None) -> dict:
    """Thermal T_ab from the image sum, with eq. 4a as printed and with the factor 1/2.

    The terms decay as n^-4, because each is a function of n beta alone at coincidence.
    The Kummer tail assumes that, and ``decay_spread`` checks it. All outputs are COMPUTED.
    """
    n_terms = N_TERMS.value if n_terms is None else n_terms
    with mp.workdps(DPS.value):
        T = mp.mpf(T)
        beta = 1 / T
        X = [mp.mpf(0)] * 4 if X is None else [mp.mpf(v) for v in X]
        split = (THERMAL_SPLIT_RATIO.value if split_ratio is None else mp.mpf(split_ratio)) * beta
        step = (THERMAL_STEP_RATIO.value if step_ratio is None else mp.mpf(step_ratio)) * beta
        out = {}
        for name, operator in (("eq4a", eq4_stress), ("with_half", stress_tensor)):
            diag = {}
            L = _sum_operator(operator, X, beta, n_terms, split, step, diag)
            src = f"imaginary-time image sum + {name} operator"
            out[name] = {
                "rho": Quantity(L[0][0], Origin.COMPUTED, src),
                "p": [Quantity(L[i][i], Origin.COMPUTED, src) for i in (1, 2, 3)],
                "trace": Quantity(-L[0][0] + L[1][1] + L[2][2] + L[3][3], Origin.COMPUTED, src),
                "max_offdiag": Quantity(
                    max(abs(L[a][c]) for a in range(4) for c in range(4) if a != c), Origin.COMPUTED, src
                ),
                "decay_spread": Quantity(diag["decay_spread"], Origin.COMPUTED, "thermal image-sum diagnostic"),
                "tail_fraction": Quantity(diag["tail_fraction"], Origin.COMPUTED, "thermal image-sum diagnostic"),
            }
    return out
