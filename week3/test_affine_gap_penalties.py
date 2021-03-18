from week3.affine_gap_penalties import align_with_affine_gap_pen, build_graph_with_affine_gap_pen
import pytest
from week2.test_global_alignment import _parse_blosum62


def build_score_dict(match, mismatch):
    score_dict = dict()
    for i in "ACGT":
        for j in "ACGT":
            score_dict[(i, j)] = match if i == j else mismatch
    return score_dict


pen0 = build_score_dict(1, -3)
test0_input = ("GTTA", "GA", (pen0, 2, 1))
pen1 = build_score_dict(1, -5)
test1_input = ("TTT", "TT", (pen1, 3, 1))
test2_input = ("GAT", "AT", (pen1, 5, 1))
test3_input = ("CCAT", "GAT", (pen1, 2, 1))
pen4 = build_score_dict(1, -2)
test4_input = ("CAGGT", "TAC", (pen4, 3, 2))
pen5 = build_score_dict(2, -3)
test5_input = ("GTTCCAGGTA", "CAGTAGTCGT", (pen5, 3, 2))
test6_input = ("AGCTAGCCTAG", "GT", (pen0, 1, 1))
pen7 = build_score_dict(2, -1)
test7_input = ("AA", "CAGTGTCAGTA", (pen7, 2, 1))
pen8 = build_score_dict(5, -2)
test8_input = ("ACGTA", "ACT", (pen8, 15, 5))


def _parse_blosum62():
    # default value when there is deletion
    blosum62 = dict()
    with open("week2/data/BLOSUM62.txt") as fd:
        lines = fd.readlines()
    w = lines[0].strip().split()
    for v_raw in lines[1:]:
        v = v_raw.strip().split()
        aa_v = v[0]
        for i, aa_w in enumerate(w):
            blosum62[(aa_w, aa_v)] = int(v[i + 1].strip())
    return blosum62


@pytest.mark.parametrize(
    "v, w, weights, direct_graph, score",
    [
        (
            *test0_input,
            [
                [None, "r", "r"],
                ["d", "dg", "r"],
                ["d", "d", "dg"],
                ["d", "d", "d"],
                ["d", "d", "dg"],
            ],
            -1,
        )
    ],
)
def test_build_graph(v, w, weights, direct_graph, score):
    d_graph, my_score = build_graph_with_affine_gap_pen(v, w, weights)
    assert my_score == score
    assert d_graph == direct_graph


@pytest.mark.parametrize(
    "v, w, weights, score, aligment",
    [
        (*test0_input, -1, ("GTTA", "G--A")),
        (*test1_input, -1, ("TTT", "TT-")),
        (*test2_input, -3, ("GAT", "-AT")),
        (*test3_input, -3, ("-CCAT", "G--AT")),
        (*test4_input, -8, ("CAGGT", "TAC--")),
        (*test5_input, -8, ("--GT--TCCAGGTA", "CAGTAGTC---GT-")),
        (*test6_input, -7, ("AGCTAGCCTAG", "-G-T-------")),
        (*test7_input, -7, ("-------A--A", "CAGTGTCAGTA")),
        (*test8_input, -12, ("ACGTA", "AC-T-")),
    ],
)
def test_align_with_affine_gap_pen(v, w, weights, score, aligment):
    m_score, m_alignment = align_with_affine_gap_pen(v, w, weights)
    assert m_score == score
    assert m_alignment == aligment


blosum62 = _parse_blosum62()
large_data = (
    "YHFDVPDCWAHRYWVENPQAIAQMEQICFNWFPSMMMKQPHVFKVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE",
    "YHEDVAHEDAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPIISATCARMRVRTVWE",
    (blosum62, 11, 1),
)


def test_align_with_affine_gap_pen_large():
    m_score, m_alignment = align_with_affine_gap_pen(*large_data)
    v = large_data[0]
    w = large_data[1]
    assert m_score == 144
    assert m_alignment[0].replace("-", "") == v
    assert m_alignment[1].replace("-", "") == w