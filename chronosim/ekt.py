"""Equations of Echeverria, Klinkhammer & Thorne (1991), as printed.

Every function here is a transcription (two passes, docs/equation_transcription.md) evaluated
exactly: the value is DERIVED from CITED formulas. Solutions of the coupled equations 4.2-4.3,
found by root finding, are COMPUTED. Nothing here is corrected to agree with the independent
model in billiards.py; disagreements are reported, not fixed.

Units of the paper (p. 1079): D = T_d = 1, B = b/D, R = r/D, speeds in D/T_d.
Angles theta, phi of the initial trajectory are those of Fig. 11; h = -B sin(theta) and
psi_A = theta + phi (p. 1089).
"""

from __future__ import annotations

import mpmath as mp

from .provenance import Origin, Quantity
from .sources import EKT_1991


def _q(value, locator, origin=Origin.DERIVED):
    return Quantity(value, origin, EKT_1991.ref(locator))


# Sec. II -------------------------------------------------------------------------------------

def traversal_psi_out(psi_in, phi):
    """Eq. 2.1b: psi_out = psi_in - 2 phi."""
    return _q(psi_in - 2 * phi, "eq. 2.1b")


def velocity_exchange_speed(B):
    """Text of p. 1083 (Fig. 6b): v1' = (D - 2B)/T_d = 1 - 2B."""
    return _q(1 - 2 * mp.mpf(B), "p. 1083, text before Fig. 6")


# Sec. III ------------------------------------------------------------------------------------

def fig7b_psi(v1):
    """Eq. 3.1: tan psi = v1 (one traversal, B << 1)."""
    return _q(mp.atan(v1), "eq. 3.1")


def fig7c_psi(v1):
    """Eq. 3.2: sin psi + tan psi = 2 v1 (two traversals, B << 1). The root in (0, pi/2) is
    unique because the left side increases from 0 to infinity there."""
    psi = mp.findroot(lambda p: mp.sin(p) + mp.tan(p) - 2 * v1, (mp.mpf(0), mp.pi / 2 - mp.mpf(10) ** -30),
                      solver="anderson")
    return _q(psi, "eq. 3.2")


def fig7c_event_times(psi):
    """Text after eq. 3.2 (p. 1084): the five event times, in order (i) to (v)."""
    c = mp.cos(psi)
    return _q([-1 / (1 + c), -c / (1 + c), mp.mpf(0), c / (1 + c), 1 / (1 + c)], "p. 1084, text after eq. 3.2")


def multitraversal_n(L, h, v1):
    """Eq. 3.7: n = 2L / (v1 - 1 - (1/2)(h/L)^2 v1)."""
    return _q(2 * L / (v1 - 1 - (h / L) ** 2 * v1 / 2), "eq. 3.7")


def multitraversal_n_min(h, v1):
    """Eq. 3.8: n_min = 3 sqrt(3 v1 / 2) |h| / (v1 - 1)^(3/2)."""
    return _q(3 * mp.sqrt(3 * v1 / 2) * abs(h) / (v1 - 1) ** mp.mpf(1.5), "eq. 3.8")


def multitraversal_speed(L, h, v1):
    """Eq. 3.6: v2 = v1' = v1 cos(h/L) (first form)."""
    return _q(v1 * mp.cos(h / L), "eq. 3.6")


# Sec. IV -------------------------------------------------------------------------------------

def initial_line(theta, phi, B):
    """p. 1089: h = -B sin(theta), psi_A = theta + phi."""
    return _q((-B * mp.sin(theta), theta + phi), "p. 1089, text")


def v1_range(theta, phi, B, R):
    """Eq. 4.8: v1min, v1max = (cos th / cos ph)(1 - 2B cos ph) -/+ 2R / cos ph."""
    base = mp.cos(theta) / mp.cos(phi) * (1 - 2 * B * mp.cos(phi))
    return _q((base - 2 * R / mp.cos(phi), base + 2 * R / mp.cos(phi)), "eq. 4.8")


def v1_from_lambda(lam, theta, phi, B, R):
    """Eq. 4.10 (and 4.11 with lam2): v = v1min + lam (v1max - v1min)."""
    lo, hi = v1_range(theta, phi, B, R).value
    return _q(lo + lam * (hi - lo), "eqs. 4.8, 4.10")


def lambda_from_speed(v, theta, phi, B, R):
    """Inverse of eqs. 4.10/4.11."""
    lo, hi = v1_range(theta, phi, B, R).value
    return _q((v - lo) / (hi - lo), "eqs. 4.8, 4.10, 4.11")


def in_dangerous_triangle(theta, phi, B):
    """Eq. 4.7, angular part: B/2 < phi < pi/2 and B - phi < theta < phi."""
    return B / 2 < phi < mp.pi / 2 and B - phi < theta < phi


def rho_of(a, b):
    """Eq. 4.4a."""
    return mp.sin((a + b) / 2) / mp.sin((a - b) / 2)


def d_of(a, b, s, R):
    """Eq. 4.4b: d = 2 s R / (1 + rho^2 + 2 rho cos a)^(1/2)."""
    rho = rho_of(a, b)
    return 2 * s * R / mp.sqrt(1 + rho ** 2 + 2 * rho * mp.cos(a))


def eq42_43(a, b, theta, phi, v1, B, R, s):
    """Residuals (left minus right) of eqs. 4.2 and 4.3 as printed."""
    th, ph = theta, phi
    d = d_of(a, b, s, R)
    h = (a - b) / 2
    r42 = (B * mp.sin(a) * mp.sin(h) * (mp.sin(th) - mp.sin(th + ph - h))
           + mp.sin(b) * mp.sin(h) * (mp.sin(th + ph) - B * (mp.sin(th) + mp.sin(th + ph + h)))
           - (v1 + d) * mp.sin((a + b) / 2) * mp.sin(a) * mp.sin(b))
    r43 = (B * mp.sin((a + b) / 2 - th - ph) * (mp.sin(a) + mp.sin(b)) + B * mp.sin(th) * mp.sin(a - b)
           - mp.sin(b) * mp.sin(a - th - ph) + d * mp.sin(a) * mp.sin(b))
    return r42, r43


def is_spurious(a, b, phi, s, R, tol):
    """Eqs. 4.6a-c: the five listed points, rho < 0, or sign(beta) != s."""
    for pa, pb in ((0, 0), (mp.pi, 0), (0, mp.pi), (mp.pi, mp.pi), (2 * phi, 0)):
        if abs(a - pa) < tol and abs(b - pb) < tol:
            return True
    if rho_of(a, b) < 0:
        return True
    return mp.sign(b) != s


def solve_eq42_43(theta, phi, v1, B, R, s, start):
    """Root of eqs. 4.2-4.3 from a starting point (a, b); COMPUTED. Returns (a, b) or None."""
    try:
        a, b = mp.findroot(lambda a_, b_: eq42_43(a_, b_, theta, phi, v1, B, R, s), start,
                           tol=mp.mpf(10) ** (-mp.mp.dps + 10))
    except (ValueError, ZeroDivisionError):
        return None
    return a, b


def central_solution(theta, phi, lam1, B, R, s):
    """Eqs. 4.12a-c (s = +1) and 4.13a-c (s = -1): (beta, eps, lam2)."""
    th, ph = theta, phi
    k = mp.tan(ph) ** 2 - mp.tan(th) ** 2
    if s == 1:
        beta = 8 * mp.sin(ph) / (mp.cos(th) * k) * B * R * lam1
        eps = 8 * mp.cos(ph) / mp.sin(th + ph) * lam1 * R
        lam2 = (1 + 2 * mp.cos(ph) / k * B) * lam1
        loc = "eqs. 4.12a-c"
    else:
        beta = -8 * mp.sin(ph) / (mp.cos(th) * k) * B * R * (1 - lam1)
        eps = -8 * mp.cos(ph) / mp.sin(th + ph) * (1 - lam1) * R
        lam2 = 1 - (1 + 2 * mp.cos(ph) / k * B) * (1 - lam1)
        loc = "eqs. 4.13a-c"
    return _q((beta, eps, lam2), loc)


def b20_first_form(theta, phi, lam1, B, s):
    """Eq. B20, first printed form: lam2 - sigma = [1 + cos^3 ph cos^2 th / (2 sin(ph-th) sin(ph+th)) B](lam1 - sigma)."""
    sigma = (1 - s) // 2
    f = 1 + mp.cos(phi) ** 3 * mp.cos(theta) ** 2 / (2 * mp.sin(phi - theta) * mp.sin(phi + theta)) * B
    return _q(sigma + f * (lam1 - sigma), "eq. B20, first form")


def b20_second_form(theta, phi, lam1, B, s):
    """Eq. B20, second printed form: lam2 - sigma = (1 + 2 cos ph / (tan^2 ph - tan^2 th) B)(lam1 - sigma)."""
    sigma = (1 - s) // 2
    f = 1 + 2 * mp.cos(phi) / (mp.tan(phi) ** 2 - mp.tan(theta) ** 2) * B
    return _q(sigma + f * (lam1 - sigma), "eq. B20, second form")


def upper_edge_414_415(phi, lam1, B, R):
    """Eq. 4.14 (beta, eps) and eq. 4.15 (beta = eps), class I, on the upper edge."""
    b414 = mp.sin(phi) * mp.cos(2 * phi) * B
    e414 = mp.cos(2 * phi) * mp.cos(phi) ** 2 / mp.sin(phi) * B
    b415 = -4 * mp.sin(phi) / mp.cos(2 * phi) * R * lam1
    return _q({"4.14": (b414, e414), "4.15": (b415, b415)}, "eqs. 4.14, 4.15")


def b30_b32(theta, phi, lam1, B, R, s=1):
    """Eqs. B21, B30, B31, B32: beta = -p/2 + sqrt(p^2/4 - q), as printed."""
    sigma = (1 - s) // 2
    mu = (phi - theta) / B
    t2 = mp.tan(phi) ** 2
    p = 2 * t2 * (2 * (lam1 - sigma) - s * mp.sin(phi)) * R - (mp.sin(phi) * mp.cos(2 * phi) - mu * t2) * B
    q = -4 * B * R * mp.sin(phi) ** 2 * (lam1 - sigma)
    disc = p ** 2 / 4 - q
    beta = -p / 2 + mp.sqrt(disc) if disc >= 0 else None
    return _q({"p": p, "q": q, "beta": beta}, "eqs. B21, B30-B32")


def appendix_b_coefficients(theta, phi, lam1, B, R, s):
    """Eqs. B3-B10 and B13-B15 as printed, for the algebraic consistency checks."""
    th, ph = theta, phi
    c, sn = mp.cos, mp.sin
    cot = lambda x: c(x) / sn(x)
    co = {
        "M": B * sn(2 * ph) * c(th) / 2,
        "N": B * (c(2 * ph) * c(th) + sn(2 * ph) * sn(th) / 4) / 2,
        "P1": -sn(ph - th) - 2 * (4 * lam1 + s - 2) * R * sn(ph)
              + B * (mp.mpf(3) / 4 * sn(2 * ph - th) - sn(2 * ph + th) / 4 - sn(th)),
        "Q1": -c(ph) * c(th) - cot(ph) * sn(ph + th) / 2 - 2 * (2 * lam1 - 1) * R * c(ph)
              + B * ((9 * c(ph) ** 2 - 1) * c(th) / 2 + sn(th) * cot(ph) * (7 * c(ph) ** 2 - 3) / 4),
        "S1": -c(2 * ph) * c(th) / c(ph) + R * (2 * (2 * lam1 - 1) / c(ph) + (4 - 8 * lam1 - s) * c(ph))
              + B * (c(th) * (mp.mpf(3) / 2 * c(ph) ** 2 - 1) + sn(2 * ph) * sn(th) / 4),
        "P2": -sn(ph - th) + 2 * s * R * sn(ph) + B * (-2 * sn(th) * c(ph) ** 2 + sn(2 * ph) * c(th) / 2),
        "Q2": -s * R * c(ph) + B * (c(th) / 2 - mp.mpf(3) / 8 * sn(2 * ph) * sn(th)),
        "S2": -c(ph - th) + s * R * c(ph) + B * (c(th) * c(ph) ** 2 + mp.mpf(5) / 4 * sn(th) * sn(2 * ph)),
        "Q": -c(ph) * c(th) - cot(ph) * sn(ph + th) / 2 - (4 * lam1 - s - 2) * R * c(ph)
             + B * (-c(th) + mp.mpf(9) / 2 * c(ph) ** 2 * c(th) + c(ph) ** 2 * cot(ph) * sn(th)),
        "S": c(th) / c(ph) - c(ph + th) + R * (2 * (1 - 2 * lam1) * c(2 * ph) / c(ph) - 2 * s * c(ph))
             + B * (-c(th) + c(th) * c(ph) ** 2 / 2 - sn(2 * ph) * sn(th)),
        "P": -4 * R * sn(ph) * (2 * lam1 + s - 1),
    }
    return _q(co, "eqs. B3-B10, B13-B15")


def b11_v1(theta, phi, lam1, B, R):
    """Eq. B11."""
    return _q((((1 - 2 * B * mp.cos(phi)) * mp.cos(theta) + 2 * R * (2 * lam1 - 1)) / mp.cos(phi)), "eq. B11")


def appendix_a_distances(a, b, theta, phi, B):
    """Eqs. A12 (theta'), A6 (a = PQ), A7 (a' = PQ'), A19 (c) and A18 (psi), as printed. Distances
    along the ball's paths, so a physical solution has a, a', c > 0 (Fig. 17)."""
    th, ph = theta, phi
    thp = ph + th - (a + b) / 2
    aL = B / mp.sin(b) * (mp.sin(th - b) - mp.sin(thp))
    aLp = aL * mp.cos((th + b + thp) / 2) / mp.cos((th - b + thp) / 2)
    c = (mp.sin(th + ph) - B * (mp.sin(th + ph + (a - b) / 2) + mp.sin(th))) / mp.sin(a)
    return _q({"theta_p": thp, "a": aL, "a_p": aLp, "c": c, "psi": -ph + (a - b) / 2}, "eqs. A6, A7, A12, A18, A19")
