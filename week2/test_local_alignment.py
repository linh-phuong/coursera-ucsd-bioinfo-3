from week2.local_alignment import build_route_with_freeride, local_alignment
import pytest
from collections import defaultdict


def build_score_dict(match, mismatch, indel):
    score_dict = defaultdict(lambda: indel)
    for i in "ACGT":
        for j in "ACGT":
            score_dict[(i, j)] = match if i == j else mismatch
    return score_dict


REPL_MAT = build_score_dict(1, -1, -2)

BT = [
    ["d", "d", "d", "d", "d"],
    ["r", "dg", "nc", "dg", "nc"],
    ["r", "nc", "dg", "r", "dg"],
    ["r", "nc", "d", "dg", "d"],
]


@pytest.mark.parametrize(
    "v, w, matrix, backtrack, maxscore, maxpos", [("GAGA", "GAT", REPL_MAT, BT, 2, (2, 2))]
)
def test_freeride_route(v, w, matrix, backtrack, maxscore, maxpos):
    bt, msc, mpo = build_route_with_freeride(v, w, matrix)
    assert bt == backtrack
    assert msc == maxscore
    assert mpo == (2, 2)


def build_score_dict(match, mismatch, indel):
    score_dict = defaultdict(lambda: indel)
    for i in "ACGT":
        for j in "ACGT":
            score_dict[(i, j)] = match if i == j else mismatch
    return score_dict


REPL_MAT1 = build_score_dict(3, -3, -1)
REPL_MAT2 = build_score_dict(1, -1, -1)
REPL_MAT3 = build_score_dict(3, -2, -1)
REPL_MAT4 = build_score_dict(2, -3, -1)


@pytest.mark.parametrize(
    "v, w, rep_score, maxscore, align",
    [
        ("GAGA", "GAT", REPL_MAT, 2, ["GA", "GA"]),
        ("AGC", "ATC", REPL_MAT1, 4, ["AG-C", "A-TC"]),
        ("AT", "AG", REPL_MAT2, 1, ["A", "A"]),
        ("TAACG", "ACGTG", REPL_MAT2, 3, ["ACG", "ACG"]),
        ("CAGAGATGGCCG", "ACG", REPL_MAT3, 6, ["CG", "CG"]),
        ("CTT", "AGCATAAAGCATT", REPL_MAT4, 5, ["C-TT", "CATT"]),
    ],
)
def test_local_alignment(v, w, rep_score, maxscore, align):
    msc, alm = local_alignment(v, w, rep_score)
    assert msc == maxscore
    assert alm == align


def _parse_pam250():
    # default value when there is deletion
    pam250 = defaultdict(lambda: -5)
    with open("week2/data/PAM250.txt") as fd:
        lines = fd.readlines()
    w = lines[0].strip().split()
    for v_raw in lines[1:]:
        v = v_raw.strip().split()
        aa_v = v[0]
        for i, aa_w in enumerate(w):
            pam250[(aa_w, aa_v)] = int(v[i + 1].strip())
    return pam250


def _parse_local_align_large():
    with open("week2/data/local_alignment.txt") as fd:
        rdat = fd.readlines()
        v = rdat[1].strip()
        w = rdat[2].strip()
        score = int(rdat[4].strip())
        align = [s.strip() for s in rdat[-2:]]
    return v, w, score, align


def test_local_alignment_large():
    pam250 = _parse_pam250()
    v, w, score, _ = _parse_local_align_large()
    myscore, myalign = local_alignment(v, w, pam250)
    assert score == myscore
    assert myalign[0].replace("-", "") in v
    assert myalign[1].replace("-", "") in w