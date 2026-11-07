"""Provenance-by-construction.

Every scientific number this package outputs is a :class:`Quantity` carrying an
:class:`Origin`. This is a purely theoretical project, so nothing here is
"measured"; the origins are:

CITED      a formula or constant transcribed from a published source
DERIVED    evaluated exactly (closed form, arbitrary precision) from CITED formulas
COMPUTED   produced by a numerical procedure of this package (image sums,
           finite differences, limits) — an approximation by construction
HEURISTIC  a value chosen by us without a source (e.g. a numerical step size)

Origins can only degrade through computation: combining quantities yields the
worst origin among the inputs, and there is no code path that promotes a
COMPUTED value to DERIVED or CITED.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum


class Origin(IntEnum):
    # Higher value = weaker claim. Combination takes the max.
    CITED = 0
    DERIVED = 1
    COMPUTED = 2
    HEURISTIC = 3


@dataclass(frozen=True)
class Quantity:
    value: object  # mpmath mpf (or float)
    origin: Origin
    source: str  # citation locator, or a description of the procedure / choice
    unit: str = ""

    def __post_init__(self):
        if not isinstance(self.origin, Origin):
            raise TypeError("origin must be an Origin")
        if not self.source or not str(self.source).strip():
            raise ValueError("every Quantity needs a non-empty source")

    def degrade(self, origin: Origin, source: str) -> "Quantity":
        """Re-label with a weaker (or equal) origin. Promotion raises."""
        if origin < self.origin:
            raise ValueError(f"cannot promote {self.origin.name} to {origin.name}")
        return Quantity(self.value, origin, source, self.unit)


def combine(value, inputs, source: str, unit: str = "", at_least: Origin = Origin.CITED) -> Quantity:
    """Build the result of a computation from its input quantities.

    The result's origin is the worst of the inputs' origins and ``at_least``
    (pass ``Origin.COMPUTED`` for any numerical approximation).
    """
    origins = [q.origin for q in inputs] + [at_least]
    return Quantity(value, max(origins), source, unit)


# Registry of every HEURISTIC constant used in a run, so it can be reported.
_HEURISTICS: dict = {}


def heuristic(name: str, value, why: str) -> Quantity:
    q = Quantity(value, Origin.HEURISTIC, why)
    _HEURISTICS[name] = q
    return q


def reported_heuristics() -> dict:
    return dict(_HEURISTICS)
