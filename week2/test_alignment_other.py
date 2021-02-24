from week2.alignment_other import (
    calculate_edit_distance,
    calculate_edit_distance_replacemat,
    hamming_distance,
)
import pytest
import sys
from week2.test_local_alignment import _parse_pam250
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


def test_edit_distance_large_with_repmat():
    v, w, ed_score = _parse_edit_distance_large()
    blosum62 = _parse_blosum62(5)
    pam250 = _parse_pam250(5)
    msc_b62 = calculate_edit_distance_replacemat(v, w, blosum62)
    msc_p250 = calculate_edit_distance_replacemat(v, w, pam250)
    assert msc_b62 == 404
    assert msc_p250 == 425