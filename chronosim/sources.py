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

CANDELAS_DEUTSCH_1977 = Source(
    key="Candelas & Deutsch (1977)",
    citation=(
        "Candelas, P. & Deutsch, D. (1977). On the vacuum stress induced by uniform acceleration or "
        "supporting the ether. Proceedings of the Royal Society of London A 354, 79-99."
    ),
    doi="10.1098/rspa.1977.0057",
    status=(
        "Volume and pages confirmed in Crossref only. PENDING: the article is closed-access and was not "
        "read, so no equation in it has been checked against the Unruh value used here."
    ),
)

KLINKHAMMER_1992 = Source(
    key="Klinkhammer (1992)",
    citation=(
        "Klinkhammer, G. (1992). Multiply connected spacetimes and closed timelike curves in "
        "semiclassical gravity. PhD thesis, California Institute of Technology."
    ),
    status=(
        "Read from the Caltech thesis repository (thesis.caltech.edu/6633). Eqs. 3.3, 3.7, 3.8 read as "
        "images, single reading."
    ),
)

VISSER_1993 = Source(
    key="Visser (1993)",
    citation=(
        "Visser, M. Hawking's chronology protection conjecture: singularity structure of the quantum "
        "stress-energy tensor. arXiv:hep-th/9303023."
    ),
    arxiv="hep-th/9303023",
    status="Read from the arXiv text. Journal version not checked.",
)

BIRRELL_DAVIES_1982 = Source(
    key="Birrell & Davies (1982)",
    citation="Birrell, N. D. & Davies, P. C. W. (1982). Quantum Fields in Curved Space. Cambridge University Press.",
    status=(
        "PENDING: the page on thermal Green functions could not be checked. The derivation used "
        "here is written out in docs/thermal_images.md and does not rely on this book."
    ),
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

EKT_1991 = Source(
    key="Echeverria, Klinkhammer & Thorne (1991)",
    citation=(
        "Echeverria, F., Klinkhammer, G. & Thorne, K. S. (1991). Billiard balls in wormhole spacetimes "
        "with closed timelike curves: Classical theory. Physical Review D 44, 1077-1099."
    ),
    doi="10.1103/PhysRevD.44.1077",
    status=(
        "Read from a scan of the published version. Equations transcribed twice (second pass blind): "
        "docs/equation_transcription.md."
    ),
)

ALL = (LI_GOTT_1998, HAWKING_1992, KIM_THORNE_1991, EKT_1991, MISNER_1967, KRW_1997, CANDELAS_DEUTSCH_1977, KLINKHAMMER_1992, VISSER_1993, BIRRELL_DAVIES_1982, HISCOCK_KONKOWSKI_1982)
