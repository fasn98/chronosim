"""The report is the source of the paper's tables, so its number formatting must not depend on the
machine. It does depend on the mpmath version: in mpmath 1.4.1 an mpf inside an f-string formats
as "3.96e-20", while mpmath 1.3.0 falls back to str(), which gives "0.0000000000000000000396" at
80 digits. The versions are pinned in requirements.txt; these tests fail under any other version.
"""

import re
from importlib.metadata import version
from pathlib import Path

import mpmath as mp
import pytest

from chronosim.config import DPS

ROOT = Path(__file__).resolve().parent.parent


def _pins():
    pins = {}
    for line in (ROOT / "requirements.txt").read_text().splitlines():
        m = re.match(r"^([A-Za-z0-9_.-]+)==(\S+)", line.strip())
        if m:
            pins[m.group(1).lower()] = m.group(2)
    return pins


@pytest.mark.parametrize("pkg", ["mpmath", "numpy", "matplotlib"])
def test_installed_versions_match_requirements(pkg):
    assert version(pkg) == _pins()[pkg], f"{pkg} {version(pkg)} installed; requirements.txt pins {_pins()[pkg]}"


@pytest.mark.parametrize("text, expected", [
    ("3.96e-20", "3.96e-20"),      # F3, finite-difference error in T^TT (outputs/report.md, main)
    ("1.64e-16", "1.64e-16"),      # F4, same quantity
])
def test_mpf_in_an_f_string_formats_like_the_published_report(text, expected):
    with mp.workdps(DPS.value):
        x = mp.mpf(text)
        assert f"{x}" == expected


def test_phase_1_and_2_report_has_no_long_fixed_point_decimals():
    """A run under another mpmath writes e.g. 0.0000000000000000000396; the committed report must not."""
    report = (ROOT / "outputs" / "report.md").read_text(encoding="utf-8")
    phases_1_2 = report.split("## Fase 3")[0]
    assert "3.96e-20" in phases_1_2 and "1.64e-16" in phases_1_2
    assert not re.search(r"\b0\.0{8,}\d", phases_1_2)
