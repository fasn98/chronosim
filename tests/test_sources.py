"""The declared state of each source must stay visible."""

from chronosim.sources import ALL, CANDELAS_DEUTSCH_1977, HISCOCK_KONKOWSKI_1982, LI_GOTT_1998


def test_li_gott_declares_the_arxiv_transcription_and_the_pending_comparison():
    assert "gr-qc/9711074" in LI_GOTT_1998.status
    assert "PENDING" in LI_GOTT_1998.status


def test_hiscock_konkowski_is_registered_as_pending():
    assert HISCOCK_KONKOWSKI_1982 in ALL
    assert HISCOCK_KONKOWSKI_1982.doi == "10.1103/PhysRevD.26.1225"
    assert "PENDING" in HISCOCK_KONKOWSKI_1982.status


def test_candelas_deutsch_declares_that_it_was_not_read():
    assert CANDELAS_DEUTSCH_1977 in ALL
    assert "PENDING" in CANDELAS_DEUTSCH_1977.status
