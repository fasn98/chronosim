"""Provenance-by-construction.

Every scientific number this package outputs is a :class:`Quantity` carrying an
:class:`Origin`. This is a purely theoretical project, so nothing here is
"measured"; the origins are:

CITED      a formula or constant transcribed from a published source
DERIVED    evaluated exactly (closed form, arbitrary precision) from CITED formulas
CERTIFIED  proven by interval arithmetic with outward rounding (mpmath.iv): a
           statement about a box, under hypotheses stated in its source, that
           holds for every point of the box. It always carries one of the
           labels of CERT_LABELS, which say what exactly was proven.
COMPUTED   produced by a numerical procedure of this package (image sums,
           finite differences, limits) — an approximation by construction
HEURISTIC  a value chosen by us without a source (e.g. a numerical step size)

Origins can only degrade through computation: combining quantities yields the
worst origin among the inputs, and there is no code path that promotes a
COMPUTED value to DERIVED or CITED. A combination never yields CERTIFIED: the
arithmetic of combine() is not part of an interval proof, so a result whose
worst input is CERTIFIED is COMPUTED. certified() is the constructor of
CERTIFIED quantities; the label rules are enforced by Quantity itself, so that
a direct construction cannot bypass them (user decision of 2026-11-07,
docs/phase3_plan.md).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum


class Origin(IntEnum):
    # Higher value = weaker claim. Combination takes the max.
    CITED = 0
    DERIVED = 1
    CERTIFIED = 2
    COMPUTED = 3
    HEURISTIC = 4


# What a CERTIFIED quantity proves; distinct labels, all approved by the user: a root certified by
# Krawczyk, always naming its system of equations (decision of 2026-11-09: "raiz certificada" alone
# is refused), a physical solution (with both interval proofs of PHYSICAL_PROOFS), and an interval
# exclusion (no root of the stated equations in the stated box; approved on 2026-11-08).
# CERT_ROOT_DIR (decision of 2026-11-10): a root of the direction machine's equations (certify.conditions_dir)
# that the independent verifier rejects; the band control tests the machine against roots of its own
# equations, not against physical solutions.
CERT_ROOT_PREFIX = "raiz certificada"
CERT_ROOT = CERT_ROOT_PREFIX + " das eqs. 4.2–4.3"
CERT_ROOT_MODEL = CERT_ROOT_PREFIX + " do modelo geométrico (billiards.py)"
CERT_ROOT_DIR = CERT_ROOT_PREFIX + " das equações da máquina por direção"
CERT_PHYSICAL = "solução física certificada"
CERT_EXCLUSION = "exclusão certificada"
CERT_LABELS = (CERT_ROOT, CERT_ROOT_MODEL, CERT_ROOT_DIR, CERT_PHYSICAL, CERT_EXCLUSION)
# CERT_PHYSICAL needs both, proven by intervals: (a) and (b) of the user decision.
PHYSICAL_PROOFS = ("positividade de a, a′ e c", "ausência de colisões extras")


@dataclass(frozen=True)
class Quantity:
    value: object  # mpmath mpf (or float)
    origin: Origin
    source: str  # citation locator, or a description of the procedure / choice
    unit: str = ""
    label: str = ""  # what was proven; set on CERTIFIED quantities only
    proofs: tuple = ()  # interval proofs behind the label (CERT_PHYSICAL needs PHYSICAL_PROOFS)

    def __post_init__(self):
        if not isinstance(self.origin, Origin):
            raise TypeError("origin must be an Origin")
        if not self.source or not str(self.source).strip():
            raise ValueError("every Quantity needs a non-empty source")
        if (self.origin is Origin.CERTIFIED) != bool(self.label):
            raise ValueError("a CERTIFIED quantity needs a label from CERT_LABELS, and only it has one")
        if self.label and self.label not in CERT_LABELS:
            raise ValueError(f"unknown certification label {self.label!r}")
        if self.label == CERT_PHYSICAL and not set(PHYSICAL_PROOFS) <= set(self.proofs):
            raise ValueError(f"{CERT_PHYSICAL!r} needs interval proofs of: {'; '.join(PHYSICAL_PROOFS)}")

    def degrade(self, origin: Origin, source: str) -> "Quantity":
        """Re-label with a weaker (or equal) origin. Promotion raises, and so does a re-label into
        CERTIFIED: only certified() makes one."""
        if origin < self.origin:
            raise ValueError(f"cannot promote {self.origin.name} to {origin.name}")
        if origin is Origin.CERTIFIED and self.origin is not Origin.CERTIFIED:
            raise ValueError("only certified() makes a CERTIFIED quantity")
        keep = origin is Origin.CERTIFIED
        return Quantity(self.value, origin, source, self.unit, self.label if keep else "", self.proofs if keep else ())


def certified(value, label: str, source: str, proofs=(), unit: str = "") -> Quantity:
    """A CERTIFIED quantity. ``source`` states the box, the equations and the hypotheses of the
    proof. CERT_PHYSICAL is refused unless ``proofs`` contains every item of PHYSICAL_PROOFS."""
    if label not in CERT_LABELS:
        raise ValueError(f"unknown certification label {label!r}")
    return Quantity(value, Origin.CERTIFIED, source, unit, label, tuple(proofs))  # the label rules are checked there


def physical_solution(value, source: str, proofs, unit: str = "") -> Quantity:
    """A physical solution: CERT_PHYSICAL when every item of PHYSICAL_PROOFS was proven by
    intervals, otherwise COMPUTED, with a source that lists what was proven and what was not."""
    proofs = tuple(proofs)
    missing = [p for p in PHYSICAL_PROOFS if p not in proofs]
    if not missing:
        return certified(value, CERT_PHYSICAL, f"{source}; provado: {'; '.join(proofs)}", proofs, unit)
    return Quantity(value, Origin.COMPUTED, f"{source}; provado: {'; '.join(proofs) or 'nada'}; não provado: {'; '.join(missing)}", unit)


def combine(value, inputs, source: str, unit: str = "", at_least: Origin = Origin.CITED) -> Quantity:
    """Build the result of a computation from its input quantities.

    The result's origin is the worst of the inputs' origins and ``at_least``
    (pass ``Origin.COMPUTED`` for any numerical approximation). A result that would be CERTIFIED
    is COMPUTED: this arithmetic is not part of an interval proof.
    """
    if at_least is Origin.CERTIFIED:
        raise ValueError("only certified() makes a CERTIFIED quantity")
    origin = max([q.origin for q in inputs] + [at_least])
    return Quantity(value, Origin.COMPUTED if origin is Origin.CERTIFIED else origin, source, unit)


# Registry of every HEURISTIC constant used in a run, so it can be reported.
_HEURISTICS: dict = {}


def heuristic(name: str, value, why: str) -> Quantity:
    q = Quantity(value, Origin.HEURISTIC, why)
    _HEURISTICS[name] = q
    return q


def reported_heuristics() -> dict:
    return dict(_HEURISTICS)
