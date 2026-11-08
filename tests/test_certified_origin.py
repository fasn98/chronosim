"""Rules of the CERTIFIED origin (user decision of 2026-11-07, docs/phase3_plan.md, second round).

- Order: CITED < DERIVED < CERTIFIED < COMPUTED < HEURISTIC.
- "raiz certificada das eqs. 4.2–4.3": when only the Krawczyk test passed.
- "solução física certificada": only when (a) a, a′, c > 0 and (b) no extra collision along the
  whole history are also proven by interval arithmetic; without (b) the physical solution stays
  COMPUTED and says what was proven.
- Combination with COMPUTED gives COMPUTED.

These tests come before any use of the label.
"""

import pytest

from chronosim.provenance import (CERT_EXCLUSION, CERT_LABELS, CERT_PHYSICAL, CERT_ROOT, PHYSICAL_PROOFS, Origin, Quantity,
                                  certified, combine, physical_solution)


def test_order_is_cited_derived_certified_computed_heuristic():
    assert [o.name for o in sorted(Origin)] == ["CITED", "DERIVED", "CERTIFIED", "COMPUTED", "HEURISTIC"]


def test_the_two_approved_labels_are_spelled_as_approved():
    assert CERT_ROOT == "raiz certificada das eqs. 4.2–4.3"
    assert CERT_PHYSICAL == "solução física certificada"
    assert set(CERT_LABELS) == {CERT_ROOT, CERT_PHYSICAL, CERT_EXCLUSION}


def test_a_certified_root_carries_its_label():
    q = certified(1, CERT_ROOT, "Krawczyk, box [..]")
    assert q.origin is Origin.CERTIFIED and q.label == CERT_ROOT


def test_certified_needs_a_known_label():
    with pytest.raises(ValueError):
        certified(1, "certificado", "x")
    with pytest.raises(ValueError):
        Quantity(1, Origin.CERTIFIED, "x")          # no label
    with pytest.raises(ValueError):
        Quantity(1, Origin.COMPUTED, "x", label=CERT_ROOT)  # a label only on CERTIFIED


@pytest.mark.parametrize("proofs", [(), (PHYSICAL_PROOFS[0],), (PHYSICAL_PROOFS[1],), ("Krawczyk",)])
def test_physical_label_needs_both_interval_proofs(proofs):
    with pytest.raises(ValueError):
        certified(1, CERT_PHYSICAL, "x", proofs=proofs)


def test_physical_label_with_both_proofs():
    q = certified(1, CERT_PHYSICAL, "x", proofs=PHYSICAL_PROOFS)
    assert q.origin is Origin.CERTIFIED and q.label == CERT_PHYSICAL


def test_physical_solution_without_the_collision_proof_stays_computed_and_says_what_was_proven():
    q = physical_solution(1, "sol", proofs=("Krawczyk", PHYSICAL_PROOFS[0]))
    assert q.origin is Origin.COMPUTED and q.label == ""
    assert "provado: Krawczyk; " + PHYSICAL_PROOFS[0] in q.source
    assert "não provado: " + PHYSICAL_PROOFS[1] in q.source


def test_physical_solution_with_every_proof_is_certified():
    q = physical_solution(1, "sol", proofs=("Krawczyk",) + PHYSICAL_PROOFS)
    assert q.origin is Origin.CERTIFIED and q.label == CERT_PHYSICAL


def test_combination_with_computed_is_computed():
    c = certified(1, CERT_ROOT, "x")
    k = Quantity(2, Origin.COMPUTED, "y")
    assert combine(3, [c, k], "z").origin is Origin.COMPUTED


def test_combination_of_certified_inputs_is_not_itself_certified():
    """Arithmetic outside the interval proof proves nothing: the result is at best COMPUTED."""
    c = certified(1, CERT_ROOT, "x")
    d = Quantity(2, Origin.CITED, "eq. 1")
    assert combine(3, [c, d], "z").origin is Origin.COMPUTED
    assert combine(3, [c, c], "z").origin is Origin.COMPUTED
    with pytest.raises(ValueError):
        combine(3, [d], "z", at_least=Origin.CERTIFIED)


def test_combination_with_heuristic_is_heuristic():
    c = certified(1, CERT_ROOT, "x")
    h = Quantity(2, Origin.HEURISTIC, "h")
    assert combine(3, [c, h], "z").origin is Origin.HEURISTIC


def test_certified_degrades_to_computed_and_loses_its_label():
    q = certified(1, CERT_ROOT, "x").degrade(Origin.COMPUTED, "y")
    assert q.origin is Origin.COMPUTED and q.label == ""


def test_nothing_is_promoted_or_relabelled_into_certified():
    with pytest.raises(ValueError):
        Quantity(1, Origin.COMPUTED, "x").degrade(Origin.CERTIFIED, "y")   # promotion
    with pytest.raises(ValueError):
        Quantity(1, Origin.DERIVED, "x").degrade(Origin.CERTIFIED, "y")    # only certified() makes one
    with pytest.raises(ValueError):
        Quantity(1, Origin.HEURISTIC, "x").degrade(Origin.CERTIFIED, "y")


def test_existing_order_relations_still_hold():
    assert Origin.CITED < Origin.DERIVED < Origin.COMPUTED < Origin.HEURISTIC
    assert combine(1, [Quantity(1, Origin.CITED, "a"), Quantity(1, Origin.DERIVED, "b")], "c").origin is Origin.DERIVED


def test_the_constructor_cannot_bypass_the_physical_proofs():
    with pytest.raises(ValueError):
        Quantity(1, Origin.CERTIFIED, "x", label=CERT_PHYSICAL)
    with pytest.raises(ValueError):
        Quantity(1, Origin.CERTIFIED, "x", label=CERT_PHYSICAL, proofs=(PHYSICAL_PROOFS[0],))
    assert Quantity(1, Origin.CERTIFIED, "x", label=CERT_PHYSICAL, proofs=PHYSICAL_PROOFS).label == CERT_PHYSICAL


def test_the_aim_and_nu_pads_are_reported_heuristics():
    from chronosim import certify
    from chronosim.provenance import reported_heuristics
    h = reported_heuristics()
    assert h["ekt_cert_aim_pad"].value > 3.1415926535897932 and h["ekt_cert_nu_pad"].value > 0
