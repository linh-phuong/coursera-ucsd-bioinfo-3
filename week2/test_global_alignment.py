from collections import defaultdict
from week2.global_alignment import (
    build_route_replace_matrix,
    build_route_with_pen,
    global_aligment,
    global_alignment_repmat,
)
import pytest


def test_build_route_with_with_pen():
    v = "ABC"
    w = "ABC"
    pen = (1, 1, 1)
    bt, sc = build_route_with_pen(v, w, pen)
    assert bt == [
        [None, "r", "r", "r"],
        ["d", "dg", "r", "r"],
        ["d", "d", "dg", "r"],
        ["d", "d", "d", "dg"],
    ]
    assert sc == 3
    v1 = "AB"
    w1 = "ABC"
    bt1, sc1 = build_route_with_pen(v1, w1, pen)
    assert bt1 == [[None, "r", "r"], ["d", "dg", "r"], ["d", "d", "dg"], ["d", "d", "d"]]
    assert sc1 == 1


@pytest.mark.parametrize(
    "v, w, pen, expmatch",
    [
        ("ABC", "ABC", (1, 1, 1), (3, ("ABC", "ABC"))),
        ("AB", "ABC", (1, 1, 1), (1, ("AB-", "ABC"))),
        ("GAGA", "GAT", (1, 1, 2), (-1, ("GAGA", "GAT-"))),
        ("GAGT", "GAT", (1, 1, 2), (1, ("GAGT", "GA-T"))),
        ("ACG", "ACT", (1, 3, 1), (0, ("ACG-", "AC-T"))),
        ("ACG", "ACT", (1, 1, 2), (1, ("ACG", "ACT"))),
        ("AT", "AG", (1, 1, 1), (0, ("AT", "AG"))),
        ("TCA", "CA", (2, 5, 1), (3, ("TCA", "-CA"))),
        ("TTTTCCTT", "CC", (1, 10, 1), (-4, ("TTTTCCTT", "----CC--"))),
        ("ACAGATTAG", "T", (2, 3, 2), (-14, ("ACAGATTAG", "-----T---"))),
        ("G", "ACATACGATG", (3, 1, 2), (-15, ("------G---", "ACATACGATG"))),
    ],
)
def test_global_aligment(v, w, pen, expmatch):
    assert expmatch == global_aligment(v, w, pen)


def _parse_blosum62(indel_punish):
    # default value when there is deletion
    blosum62 = defaultdict(lambda: -indel_punish)
    with open("week2/data/BLOSUM62.txt") as fd:
        lines = fd.readlines()
    w = lines[0].strip().split()
    for v_raw in lines[1:]:
        v = v_raw.strip().split()
        aa_v = v[0]
        for i, aa_w in enumerate(w):
            blosum62[(aa_w, aa_v)] = int(v[i + 1].strip())
    return blosum62


def _parse_global_align_large():
    with open("week2/data/global_alignment.txt") as fd:
        rdat = fd.readlines()
        v = rdat[1].strip()
        w = rdat[2].strip()
        score = int(rdat[4].strip())
        align = [s.strip() for s in rdat[-2:]]
    return v, w, score, align


BLOSUM62 = _parse_blosum62(5)
DM = [
    [None, "r", "r", "r", "r"],
    ["d", "dg", "r", "r", "r"],
    ["d", "d", "dg", "r", "r"],
    ["d", "d", "d", "dg", "r"],
    ["d", "d", "d", "d", "dg"],
]
DM1 = [
    [None, "r", "r", "r", "r"],
    ["d", "dg", "r", "r", "r"],
    ["d", "d", "dg", "r", "r"],
    ["d", "d", "d", "dg", "dg"],
]


@pytest.mark.parametrize(
    "v, w, matrix, dmat, score",
    [("ILYP", "ILIP", BLOSUM62, DM, 14), ("ILYP", "ILP", BLOSUM62, DM1, 10)],
)
def test_build_route_replace_matrix(v, w, matrix, dmat, score):
    assert build_route_replace_matrix(v, w, matrix) == (dmat, score)


@pytest.mark.parametrize(
    "v, w, matrix, align, score",
    [
        ("ILYP", "ILIP", BLOSUM62, ("ILYP", "ILIP"), 14),
        ("ILYP", "ILP", BLOSUM62, ("ILYP", "IL-P"), 10),
    ],
)
def test_global_align_blosum(v, w, matrix, align, score):
    assert global_alignment_repmat(v, w, matrix) == (score, align)


def test_global_align_large():
    v, w, score, align = _parse_global_align_large()
    myscore, myalign = global_alignment_repmat(v, w, BLOSUM62)
    assert myscore == score
    assert myalign[0] == align[0]
    assert myalign[1] == align[1]
