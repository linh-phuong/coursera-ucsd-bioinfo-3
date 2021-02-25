from week2.alignment_other import (
    build_route_with_freeride_to_first_match,
    calculate_edit_distance,
    calculate_edit_distance_replacemat,
    find_fitting_alignment,
    find_fitting_positions,
    find_overlap_alignment,
    find_overlap_positions,
    hamming_distance,
)
import pytest
import sys
from week2.test_local_alignment import _parse_pam250, build_score_dict
from week2.test_global_alignment import _parse_blosum62

sys.setrecursionlimit(1500)


@pytest.mark.parametrize("v, w, distance", [("ABC", "AB-", 1), ("GAGA", "GAT-", 2)])
def test_hamming_distance(v, w, distance):
    assert hamming_distance(v, w) == distance


@pytest.mark.parametrize(
    "v, w, pen, edit_distance",
    [
        ("ABC", "AB", (1, 1, 1), 1),
        ("ABC", "AB", (0, 1, 1), 1),
        ("GAGA", "GAT", (1, 1, 1), 2),
        ("GAGA", "GAT", (0, 1, 1), 2),
        ("AC", "AC", (1, 1, 1), 0),
        ("AT", "G", (1, 1, 1), 2),
        ("CAGACCGAGTTAG", "CGG", (0, 1, 1), 10),
        ("CGT", "CAGACGGTGACG", (1, 1, 1), 9),
    ],
)
def test_edit_distance(v, w, pen, edit_distance):
    assert edit_distance == calculate_edit_distance(v, w, pen)


def _parse_edit_distance_large():
    with open("week2/data/edit_distance.txt") as fd:
        fd.readline()
        v = fd.readline().strip()
        w = fd.readline().strip()
        fd.readline()
        ed_score = int(fd.readline().strip())
    return v, w, ed_score


def test_edit_distance_large():
    v, w, ed_score = _parse_edit_distance_large()
    msc = calculate_edit_distance(v, w, (0, 1, 1))
    assert msc == ed_score


# This test shows that, using the blosum and the pam matrix
# we do not obtain the alignment score in the test data
def test_edit_distance_large_with_repmat():
    v, w, _ = _parse_edit_distance_large()
    blosum62 = _parse_blosum62(5)
    pam250 = _parse_pam250(5)
    msc_b62 = calculate_edit_distance_replacemat(v, w, blosum62)
    msc_p250 = calculate_edit_distance_replacemat(v, w, pam250)
    assert msc_b62 == 404
    assert msc_p250 == 425


WEIGHT = build_score_dict(1, -1, -1)
WEIGHT1 = build_score_dict(1, -1, -2)
WEIGHT2 = build_score_dict(1, -5, -1)
WEIGHT3 = build_score_dict(2, -3, -1)
WEIGHT4 = build_score_dict(10, -1, -1)
WEIGHT5 = build_score_dict(1, -2, -2)
WEIGHT6 = build_score_dict(1, -1, -5)
WEIGHT7 = build_score_dict(3, -2, -1)

BT = [
    ["fr", "fr", "fr", "fr", "fr", "fr"],
    ["d", "dg", "r", "d", "d", "dg"],
    ["d", "d", "dg", "r", "r", "d"],
    ["d", "d", "d", "dg", "r", "dg"],
]
AL = [[0, 0, 0, 0, 0, 0], [-1, 1, 0, -1, -1, 1], [-2, 0, 2, 1, 0, 0], [-3, -1, 1, 1, 0, 1]]

BT1 = [
    ["fr", "fr", "fr", "fr", "fr", "fr"],
    ["d", "dg", "r", "d", "d", "dg"],
    ["d", "d", "d", "dg", "r", "d"],
    ["d", "d", "d", "d", "d", "dg"],
]
AL1 = [[0, 0, 0, 0, 0, 0], [-1, 1, 0, -1, -1, 1], [-2, 0, -1, 1, 0, 0], [-3, -1, -2, 0, -1, 1]]

BT2 = [["fr", "fr", "fr"], ["d", "d", "dg"], ["d", "d", "d"]]
AL2 = [[0, 0, 0], [-2, -2, 1], [-4, -4, -1]]

BT3 = [
    ["fr", "fr", "fr", "fr"],
    ["d", "d", "dg", "r"],
    ["d", "d", "d", "dg"],
    ["d", "dg", "d", "d"],
]
AL3 = [[0, 0, 0, 0], [-2, -2, 1, -1], [-4, -4, -1, 2], [-6, -3, -3, 0]]


@pytest.mark.parametrize(
    "v, w, weight, backtrack, align",
    [
        ("TAGGT", "TAT", WEIGHT, BT, AL),
        ("TAGGT", "TGT", WEIGHT2, BT1, AL1),
        ("GA", "AT", WEIGHT5, BT2, AL2),
        ("GAT", "ATG", WEIGHT5, BT3, AL3),
    ],
)
def test_build_route_with_freeride_to_first_match(v, w, weight, backtrack, align):
    bt, sc = build_route_with_freeride_to_first_match(v, w, weight)
    assert bt == backtrack
    assert sc == align


@pytest.mark.parametrize(
    "align, score, positions",
    [(AL, 1, [(3, 3), (3, 5)])],
)
def test_find_fitting_positions(align, score, positions):
    sc, pos = find_fitting_positions(align)
    assert sc == score
    assert pos == positions


# There are many possible solution to this function
# The test check the score of the alignment and whether
# the aligment when remove "-" is exactly w
# but here I also want to show the possible alignments
# hence there is the input alignment result which is my abitrary result
@pytest.mark.parametrize(
    "v, w, weight, score, alignment",
    [
        ("TAGGT", "TAT", WEIGHT, 1, [("TAG", "TAT"), ("TAGGT", "TA--T")]),
        ("GAGA", "GAT", WEIGHT1, 1, [("GAG", "GAT")]),
        ("CCAT", "AT", WEIGHT, 2, [("AT", "AT")]),
        ("CACGTC", "AT", WEIGHT2, 0, [("A-", "AT"), ("-T", "AT")]),
        ("ATCC", "AT", WEIGHT, 2, [("AT", "AT")]),
        ("ACGACAGAG", "CGAGAGGTT", WEIGHT3, 7, [("CGACAGAG---", "CGA--GAGGTT")]),
        ("CAAGACTACTATTAG", "GG", WEIGHT4, 10, [("GACTACTATTAG", "G----------G")]),
    ],
)
def test_find_fitting_alignment(v, w, weight, score, alignment):
    sc, al_str = find_fitting_alignment(v, w, weight)
    assert sc == score
    assert al_str == alignment
    for al in al_str:
        w_sb = al[1]
        assert w_sb.replace("-", "") == w


def _parse_fitting_alignment_large():
    with open("week2/data/fitting_alignment.txt") as fd:
        fd.readline()
        v = fd.readline().strip()
        w = fd.readline().strip()
        fd.readline()
        score = int(fd.readline().strip())
        v_al = fd.readline().strip()
        w_al = fd.readline().strip()
    return v, w, score, v_al, w_al


def test_fitting_alignment_large():
    v, w, score, v_al, w_al = _parse_fitting_alignment_large()
    m_sc, m_al = find_fitting_alignment(v, w, WEIGHT)
    assert m_sc == score
    assert m_al[0][0] == v_al
    assert m_al[0][1] == w_al


@pytest.mark.parametrize("align, score, positions", [(AL2, 1, [(1, 2)]), (AL3, 2, [(2, 3)])])
def test_find_overlap_positions(align, score, positions):
    sc, pos = find_overlap_positions(align)
    assert sc == score
    assert pos == positions


@pytest.mark.parametrize(
    "v, w, weight, score, alignments",
    [
        ("GA", "AT", WEIGHT5, 1, [("A", "A")]),
        ("GAT", "ATG", WEIGHT5, 2, [("AT", "AT")]),
        ("GAGA", "GAT", WEIGHT1, 2, [("GA", "GA")]),
        ("CCAT", "AT", WEIGHT, 2, [("AT", "AT")]),
        ("GAT", "CAT", WEIGHT2, 1, [("-AT", "CAT")]),
        ("ATCACT", "AT", WEIGHT2, 1, [("ACT", "A-T")]),
        ("ATCACT", "ATG", WEIGHT6, 0, [("CT", "AT")]),
        ("CAGAGATGGCCG", "ACG", WEIGHT7, 5, [("ATGGCCG", "A---C-G")]),
        ("CTT", "AGCATAAAGCATT", WEIGHT3, 0, [("--C-TT", "AGCAT-")]),
    ],
)
def test_find_overlap_alignment(v, w, weight, score, alignments):
    sc, al = find_overlap_alignment(v, w, weight)
    assert sc == score
    assert al == alignments


def _parse_overlap_alignment_large():
    with open("week2/data/overlap_alignment.txt") as fd:
        fd.readline()
        v = fd.readline().strip()
        w = fd.readline().strip()
        fd.readline()
        score = int(fd.readline().strip())
        v_suff = fd.readline().strip()
        w_pref = fd.readline().strip()
    return v, w, score, v_suff, w_pref


# The large test uses WEIGHT5
def test_overlap_alignment_large():
    v, w, score, v_suff, w_pref = _parse_overlap_alignment_large()
    msc, m_align = find_overlap_alignment(v, w, WEIGHT5)
    assert msc == score
    mv_suff = m_align[0][0].replace("-", "")
    mw_pref = m_align[0][1].replace("-", "")
    len_suff = len(mv_suff)
    len_pref = len(mw_pref)
    assert mv_suff == v[-len_suff:]
    assert len_suff == len(v_suff.replace("-", ""))
    assert mw_pref == w[0:len_pref]
    assert len_pref == len(w_pref.replace("-", ""))