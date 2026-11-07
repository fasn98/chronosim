"""Numerical settings shared by all modules (all HEURISTIC, reported on every run)."""

from .provenance import heuristic

DPS = heuristic("dps", 80, "working precision in decimal digits (mpmath)")
