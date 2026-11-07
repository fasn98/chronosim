"""F7: arithmetic check of the order-of-magnitude estimates of Kim & Thorne (1991), Sec. II.

These are estimates with "~". Evaluating them checks only the arithmetic between the
printed inputs and the printed numbers, not whether the estimates are right. Every
output is DERIVED from CITED inputs and is labelled as an order of magnitude.
"""

from __future__ import annotations

import mpmath as mp

from .provenance import Origin, Quantity
from .sources import HAWKING_1992, KIM_THORNE_1991

OOM = "order of magnitude"

# Inputs printed in the articles (CITED).
PLANCK_LENGTH_CM = Quantity(mp.mpf("1e-33"), Origin.CITED, KIM_THORNE_1991.ref("p. 3930: l_P ~ 10^-33 cm, not checked visually; implied by l_P/D ~ 10^-35 with D = 1 m, p. 3935"), "cm")
PLANCK_MASS_G = Quantity(mp.mpf("1e-5"), Origin.CITED, KIM_THORNE_1991.ref("p. 3943: m_P ~ 10^-5 g"), "g")
SEPARATION_CM = Quantity(mp.mpf(100), Origin.CITED, KIM_THORNE_1991.ref("p. 3935: D ~ 1 m"), "cm")

# Numbers printed in the articles, to compare with.
PRINTED_LP_OVER_D = Quantity(mp.mpf("1e-35"), Origin.CITED, KIM_THORNE_1991.ref("p. 3935, after eq. 14"))
PRINTED_MAX_STRAIN = Quantity(mp.mpf("1e-35"), Origin.CITED, KIM_THORNE_1991.ref("p. 3935, after eq. 13"))
PRINTED_MAX_DENSITY = Quantity(mp.mpf("1e49"), Origin.CITED, KIM_THORNE_1991.ref("p. 3935, after eq. 14"), "g/cm^3")
PRINTED_HAWKING_DT = Quantity(mp.mpf("1e-35"), Origin.CITED, KIM_THORNE_1991.ref("p. 3937: Delta t ~ 10^-35 l_P"), "l_P")
PRINTED_PLANCK_DENSITY_HAWKING = Quantity(mp.mpf("1e94"), Origin.CITED, HAWKING_1992.ref("p. 603: 10^94 g/cc"), "g/cm^3")


def metric_perturbation(b_over_D, power, lp_over_D, lp_over_dt) -> Quantity:
    """delta g^VP ~ (b/D)^(N or N-1) (l_P/D) (l_P/Delta t), eq. 9c."""
    return Quantity(b_over_D ** power * lp_over_D * lp_over_dt, Origin.DERIVED, f"{KIM_THORNE_1991.ref('eq. 9c')}, {OOM}")


def quantum_gravity_ratio(b_over_D, power, D_over_lp) -> Quantity:
    """delta g^QG / delta g^VP ~ (D / l_P) (D / b)^(N or N-1), eq. 11."""
    return Quantity(D_over_lp / b_over_D ** power, Origin.DERIVED, f"{KIM_THORNE_1991.ref('eq. 11')}, {OOM}")


def lp_over_D() -> Quantity:
    return Quantity(PLANCK_LENGTH_CM.value / SEPARATION_CM.value, Origin.DERIVED, f"l_P / D with D = 1 m, {OOM}")


def max_strain(b_over_D=1, power=0) -> Quantity:
    """(delta L / L)_max ~ (b/D)^(N or N-1) (l_P/D), eq. 13; N - 1 = 0 near a mouth."""
    return Quantity(b_over_D ** power * lp_over_D().value, Origin.DERIVED, f"{KIM_THORNE_1991.ref('eq. 13')}, {OOM}")


def planck_density() -> Quantity:
    """m_P / l_P^3 in g/cm^3, from the inputs printed by Kim & Thorne."""
    return Quantity(PLANCK_MASS_G.value / PLANCK_LENGTH_CM.value ** 3, Origin.DERIVED, f"m_P / l_P^3, {OOM}", "g/cm^3")


def max_density() -> Quantity:
    """T^00_max ~ (l_P / D)(m_P / l_P^3), eq. 14, in g/cm^3."""
    return Quantity(lp_over_D().value * planck_density().value, Origin.DERIVED, f"{KIM_THORNE_1991.ref('eq. 14')}, {OOM}", "g/cm^3")


def hawking_cutoff_dt() -> Quantity:
    """Delta t, in units of l_P, where eq. 9c reaches 1 with N - 1 = 0: (l_P/D)(l_P/Delta t) = 1, so Delta t / l_P = l_P / D."""
    return Quantity(lp_over_D().value, Origin.DERIVED,
                    f"{KIM_THORNE_1991.ref('eq. 9c')} = 1 solved in closed form (Hawking's cutoff, Sec. II F), {OOM}", "l_P")
