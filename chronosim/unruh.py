"""Arbiter independent of Li & Gott and Kim & Thorne: the Rindler vacuum from the thermal density.

A conformal scalar at temperature T has energy density pi^2 T^4 / 30. With the local
Unruh temperature T = 1 / (2 pi xi), the Rindler vacuum is its opposite,
rho = -1 / (480 pi^2 xi^4). The relation is the one attributed to Candelas & Deutsch
(1977), whose text has not been read (see its status in ``sources.py``).
"""

from __future__ import annotations

import mpmath as mp

from .config import DPS
from .provenance import Origin, Quantity
from .sources import CANDELAS_DEUTSCH_1977


def unruh_temperature(xi):
    """T = 1 / (2 pi xi)."""
    return 1 / (2 * mp.pi * xi)


def thermal_density_closed(T):
    """pi^2 T^4 / 30."""
    return mp.pi ** 2 * T ** 4 / 30


def thermal_density_integral(T) -> Quantity:
    """Bose integral for one massless scalar degree of freedom, (1 / 2 pi^2) int k^3 / (e^(k/T) - 1) dk."""
    with mp.workdps(DPS.value):
        value = mp.quad(lambda k: k ** 3 / mp.expm1(k / T), [0, T, 10 * T, mp.inf]) / (2 * mp.pi ** 2)
    return Quantity(value, Origin.COMPUTED, "numerical Bose integral")


def rindler_density(xi) -> Quantity:
    """rho = -pi^2 T^4 / 30 at the local Unruh temperature.

    HEURISTIC: the relation "Rindler vacuum = minus thermal" is cited from Candelas & Deutsch
    (1977) without the article having been read. The thermal arbiter is the main one.
    """
    with mp.workdps(DPS.value):
        value = -thermal_density_closed(unruh_temperature(mp.mpf(xi)))
    return Quantity(value, Origin.HEURISTIC,
                    "relação vácuo de Rindler = −térmico citada de " + CANDELAS_DEUTSCH_1977.ref("artigo não lido"))
