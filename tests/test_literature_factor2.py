"""Literature check of the factor 2: Klinkhammer (1992) against Kim & Thorne (1991)."""

import mpmath as mp

from chronosim import kim_thorne as kt
from chronosim import thermal
from chronosim import wormhole as wh
from chronosim.pointsplit import stress_tensor

METRIC = (-1, 1, 1, 1)


def _klinkhammer_eq33(F, X, delta, h):
    """Klinkhammer eq. 3.3 in flat space (R = 0), as printed, by central differences, split along x."""

    def shift(P, i, s):
        P = list(P)
        P[i] += s
        return P

    def at(Xp):
        mixed = [[(F(shift(X, a, h), shift(Xp, c, h)) - F(shift(X, a, h), shift(Xp, c, -h))
                   - F(shift(X, a, -h), shift(Xp, c, h)) + F(shift(X, a, -h), shift(Xp, c, -h))) / (4 * h * h)
                  for c in range(4)] for a in range(4)]
        same = [[(F(shift(shift(X, a, h), c, h), Xp) - F(shift(shift(X, a, h), c, -h), Xp)
                  - F(shift(shift(X, a, -h), c, h), Xp) + F(shift(shift(X, a, -h), c, -h), Xp)) / (4 * h * h)
                 for c in range(4)] for a in range(4)]
        other = [[(F(X, shift(shift(Xp, a, h), c, h)) - F(X, shift(shift(Xp, a, h), c, -h))
                   - F(X, shift(shift(Xp, a, -h), c, h)) + F(X, shift(shift(Xp, a, -h), c, -h))) / (4 * h * h)
                  for c in range(4)] for a in range(4)]
        box = sum(METRIC[c] * same[c][c] for c in range(4))
        box_p = sum(METRIC[c] * other[c][c] for c in range(4))
        cross = sum(METRIC[c] * mixed[c][c] for c in range(4))
        return [[mp.mpf(1) / 6 * (mixed[c][a] + mixed[a][c]) - mp.mpf(1) / 12 * METRIC[a] * (a == c) * cross
                 - mp.mpf(1) / 12 * (same[a][c] + other[a][c]) + mp.mpf(1) / 48 * METRIC[a] * (a == c) * (box + box_p)
                 for c in range(4)] for a in range(4)]

    plus, minus = at(shift(X, 1, delta)), at(shift(X, 1, -delta))
    return [[(plus[a][c] + minus[a][c]) / 2 for c in range(4)] for a in range(4)]


def test_klinkhammer_operator_equals_the_operator_with_one_half_on_the_thermal_state():
    # One image pair of the thermal Hadamard function: symmetric in x and x', a solution of the wave equation.
    mp.mp.dps = 80
    F = lambda P, Q: thermal.hadamard_image_pair(1, P, Q, mp.mpf(1))
    X = [mp.mpf(0)] * 4
    K = _klinkhammer_eq33(F, X, mp.mpf("1e-10"), mp.mpf("1e-11"))
    H = stress_tensor(F, X)
    for a in range(4):
        assert abs(K[a][a] / H[a][a] - 1) < mp.mpf("1e-18")  # measured at most 4.0e-20
    for a in range(4):
        for c in range(4):
            if a != c:
                assert abs(K[a][c]) < mp.mpf("1e-30") * abs(H[0][0])  # measured 0


def test_klinkhammer_eq38_is_half_of_kim_thorne_eq67_and_matches_the_operator_with_one_half():
    mp.mp.dps = 80
    D, b, X, delta = mp.mpf(1), mp.mpf("0.01"), mp.mpf("0.25"), mp.mpf("1e-6")
    xi = kt.doppler_xi(mp.mpf("0.5"))
    num = wh.stress_between_numeric(X - delta, X, mp.mpf(0), D, b, "0.5", delta)
    n = num["n_terms"]
    klinkhammer = kt.klinkhammer_kk_on_axis(delta, X, D, b, xi, n).value
    kim_thorne = kt.stress_between_kk(X - delta, X, mp.mpf(0), D, b, xi, n).value
    assert abs(klinkhammer / kim_thorne - mp.mpf(1) / 2) < mp.mpf("1e-70")  # an algebraic identity
    assert abs(num["with_half"]["TT"].value / klinkhammer - 1) < mp.mpf("1e-7")  # the delta / X correction, measured 3.3e-8
