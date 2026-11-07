"""Bibliographic sources. Every CITED quantity or formula must point to one of these."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Source:
    key: str
    citation: str
    doi: str | None = None
    arxiv: str | None = None
    status: str | None = None  # what has and has not been checked against this source

    def ref(self, locator: str) -> str:
        """A precise pointer, e.g. ``LI_GOTT_1998.ref("eq. 12")``."""
        return f"{self.key}, {locator}"


LI_GOTT_1998 = Source(
    key="Li & Gott (1998)",
    citation=(
        "Li, L.-X. & Gott, J. R. (1998). Self-consistent vacuum for Misner space and the "
        "chronology protection conjecture. Physical Review Letters 80, 2980-2983."
    ),
    doi="10.1103/PhysRevLett.80.2980",
    arxiv="gr-qc/9711074",
    status=(
        "Transcribed from the arXiv version gr-qc/9711074. Comparison with the published version "
        "(PRL 80, 2980) is PENDING: no access to it."
    ),
)

HAWKING_1992 = Source(
    key="Hawking (1992)",
    citation="Hawking, S. W. (1992). Chronology protection conjecture. Physical Review D 46, 603-611.",
    doi="10.1103/PhysRevD.46.603",
    status="Read from a PDF whose title, volume and first page were checked in the PDF itself.",
)

KIM_THORNE_1991 = Source(
    key="Kim & Thorne (1991)",
    citation=(
        "Kim, S.-W. & Thorne, K. S. (1991). Do vacuum fluctuations prevent the creation of closed "
        "timelike curves? Physical Review D 43, 3929-3947."
    ),
    doi="10.1103/PhysRevD.43.3929",
    status="Read from a scan of the published version. Transcriptions: docs/equation_transcription.md.",
)

MISNER_1967 = Source(
    key="Misner (1967)",
    citation=(
        "Misner, C. W. (1967). Taub-NUT space as a counterexample to almost anything. "
        "In: Relativity Theory and Astrophysics I, Lectures in Applied Mathematics 8, 160-169."
    ),
)

KRW_1997 = Source(
    key="Kay, Radzikowski & Wald (1997)",
    citation=(
        "Kay, B. S., Radzikowski, M. J. & Wald, R. M. (1997). Quantum field theory on spacetimes "
        "with a compactly generated Cauchy horizon. Communications in Mathematical Physics 183, 533-556."
    ),
    doi="10.1007/s002200050042",
    arxiv="gr-qc/9603012",
)

HISCOCK_KONKOWSKI_1982 = Source(
    key="Hiscock & Konkowski (1982)",
    citation=(
        "Hiscock, W. A. & Konkowski, D. A. (1982). Quantum vacuum energy in Taub-NUT "
        "(Newman-Unti-Tamburino)-type cosmologies. Physical Review D 26, 1225."
    ),
    doi="10.1103/PhysRevD.26.1225",
    status=(
        "PENDING: no open version of the article is available. The phase that depends on it is blocked "
        "until the PDF is obtained. No formula from it may be used or reconstructed from memory."
    ),
)

ALL = (LI_GOTT_1998, HAWKING_1992, KIM_THORNE_1991, MISNER_1967, KRW_1997, HISCOCK_KONKOWSKI_1982)
