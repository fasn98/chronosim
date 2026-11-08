"""The CSVs that carry the phase 2 tables otherwise found only in report.md (read by
paper/build_tables.py). They must exist, say their origin, and agree with the report."""

import csv
from pathlib import Path

import pytest

OUT = Path(__file__).resolve().parent.parent / "outputs"
TABLES = ("kt_arbiter_eq4a_vs_eq11.csv", "kt_thermal_arbiter.csv", "kt_thermal_modes.csv", "kt_unruh_arbiter.csv",
          "kim_thorne_f2.csv", "kim_thorne_f5.csv", "kim_thorne_f6.csv")


def _rows(name):
    with open(OUT / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


@pytest.mark.parametrize("name", TABLES)
def test_table_exists_and_every_row_states_an_origin(name):
    rows = _rows(name)
    assert rows and all(r["origin"].startswith("COMPUTED") for r in rows)


def test_unruh_table_says_its_reference_is_heuristic():
    assert all("HEURISTIC" in r["origin"] for r in _rows("kt_unruh_arbiter.csv"))


def test_f7_keeps_printed_and_recomputed_origins_apart():
    rows = _rows("kim_thorne_f7.csv")
    assert {r["printed_origin"] for r in rows} == {"CITED"}
    assert {r["recomputed_origin"] for r in rows} <= {"DERIVED", "COMPUTED"}


def test_f6_labels_match_the_report():
    """Delta T is in Planck lengths and b = 10 l_P, so Delta T / b is a tenth of Delta T (report.md, F6)."""
    labels = [r["parameter"] for r in _rows("kim_thorne_f6.csv") if r["quantity"] == "eq54_dev"]
    report = (OUT / "report.md").read_text(encoding="utf-8")
    for label in labels:
        value = label.split("= ")[1]
        assert f"ΔT/b = {value}:" in report, label


def test_summary_has_the_scalar_results():
    keys = {r["key"] for r in _rows("summary.csv")}
    assert {"phase1_worst_rel_error", "phase1_horizon_slope", "f1_eq31", "f3_slope_eq4a", "f4_slope_eq4a", "f5_eq43_max_dev"} <= keys
