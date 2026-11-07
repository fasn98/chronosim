import pytest

from chronosim.provenance import Origin, Quantity, combine, heuristic, reported_heuristics


def test_origin_cannot_be_promoted():
    q = Quantity(1.0, Origin.COMPUTED, "image sum")
    with pytest.raises(ValueError):
        q.degrade(Origin.DERIVED, "pretend it is exact")
    assert q.degrade(Origin.HEURISTIC, "worse is allowed").origin is Origin.HEURISTIC


def test_combination_takes_the_worst_origin():
    a = Quantity(1, Origin.CITED, "paper, eq. 1")
    b = Quantity(2, Origin.COMPUTED, "numerics")
    assert combine(3, [a, b], "sum").origin is Origin.COMPUTED
    assert combine(3, [a], "x2").origin is Origin.CITED
    assert combine(3, [a], "approx", at_least=Origin.COMPUTED).origin is Origin.COMPUTED


def test_every_quantity_needs_a_source():
    with pytest.raises(ValueError):
        Quantity(1, Origin.CITED, "  ")
    with pytest.raises(TypeError):
        Quantity(1, "CITED", "paper")


def test_heuristics_are_registered_for_reporting():
    heuristic("test_knob", 3, "used only in this test")
    h = reported_heuristics()
    assert h["test_knob"].origin is Origin.HEURISTIC
    # the package's own numerical knobs are always reported
    import chronosim.pointsplit  # noqa: F401

    assert {"dps", "delta", "h", "n_terms"} <= set(reported_heuristics())
