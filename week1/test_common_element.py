from week1.common_element import (
    convert_to_backtrack,
    find_highest_score_track,
    find_longest_common,
    find_longest_route_value_MT,
    find_longest_route_value_MTD,
    find_min_elements,
    max_with_start_end,
    max_with_tie,
    min_nb_elements,
    construct_route,
    select_route_with_startend,
    track_longest_route,
)
import numpy as np
from pathlib import Path
import pytest
import textwrap
import sys
import re

sys.setrecursionlimit(1500)


def test_find_longest_common():
    assert find_longest_common("ACTGCA", "CATCGT") == "ACT"
    assert find_longest_common("ACTGCA", "CATCGC") == "ACGC"
    assert find_longest_common("AACTG", "CAAG") == "AAG"
    assert find_longest_common("ATGTTATA", "ATCGTCC") == "ATGT"


def test_min_nb_elements():
    c = [1, 4, 5]
    assert min_nb_elements(0, c) == 0
    assert min_nb_elements(1, c) == 1
    assert min_nb_elements(2, c) == 2
    assert min_nb_elements(3, c) == 3
    assert min_nb_elements(4, c) == 1
    assert min_nb_elements(22, c) == 5
    c = [50, 25, 20, 10, 5, 1]
    assert min_nb_elements(40, c) == 2


def test_find_min_elements():
    c = [1, 4, 5]
    assert find_min_elements(0, c) == []
    assert find_min_elements(1, c) == [1]
    assert find_min_elements(2, c) == [1, 1]
    assert find_min_elements(4, c) == [4]
    assert find_min_elements(22, c) == [5, 5, 4, 4, 4]


def test_find_longest_route():
    D = np.array([[1, 1, 1], [1, 1, 1]])
    R = np.array([[1, 1], [1, 1], [1, 1]])
    assert find_longest_route_value_MT(2, 2, D, R) == 4
    D = np.array([[1, 1, 1], [1, 1, 1]])
    R = np.array([[4, 1], [1, 1], [1, 1]])
    assert find_longest_route_value_MT(2, 2, D, R) == 7
    D = np.array([[1, 0, 2, 4, 3], [4, 6, 5, 2, 1], [4, 4, 5, 2, 1], [5, 6, 8, 5, 3]])
    R = np.array([[3, 2, 4, 0], [3, 2, 4, 2], [0, 7, 3, 3], [3, 3, 0, 2], [1, 3, 2, 2]])
    assert find_longest_route_value_MT(4, 4, D, R) == 34


def _parse_longest_route_test(file):
    with open(file) as fd:
        L = fd.readline().strip().split()
        n = int(L[0].strip())
        m = int(L[1].strip())
        D = []
        R = []
        while True:
            L = fd.readline().strip().split()
            if L[0] == "-":
                break
            D.append([int(i.strip()) for i in L])
        while True:
            L = fd.readline().strip().split()
            if not L:
                break
            R.append([int(i.strip()) for i in L])
        return n, m, np.array(D), np.array(R)


def _parse_longest_route_large():
    with open("week1/data/Manhattan_tourist.txt") as fd:
        fd.readline()
        L = fd.readline().strip().split()
        n = int(L[0].strip())
        m = int(L[1].strip())
        D = []
        R = []
        while True:
            L = fd.readline().strip().split()
            if L[0] == "-":
                break
            D.append([int(i.strip()) for i in L])
        while True:
            L = fd.readline().strip().split()
            if not L:
                break
            R.append([int(i.strip()) for i in L])
        fd.readline()
        r_exp = int(fd.readline().strip())
        return n, m, np.array(D), np.array(R), r_exp


def test_parse_longest_route(tmp_path):
    fn = tmp_path / "a.txt"
    fn.write_text(
        textwrap.dedent(
            """
    3 3
    0 1 2
    3 4 5
    -
    6 7
    8 9
    10 11
    """
        ).strip()
    )
    n, m, D, R = _parse_longest_route_test(fn)
    assert n == 3
    assert m == 3
    assert (D == np.array([[0, 1, 2], [3, 4, 5]])).all()
    assert (R == np.array([[6, 7], [8, 9], [10, 11]])).all()


TEST_DIR = Path("week1/data/LongestPathGrid/inputs/")


@pytest.mark.parametrize("inp", sorted(TEST_DIR.glob("*")))
def test_longest_route(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    n, m, D, R = _parse_longest_route_test(inp)
    with open(out) as fd:
        r_exp = int(fd.readline().strip())
    r = find_longest_route_value_MT(n, m, D, R)
    assert r == r_exp


def test_longest_route_large():
    n, m, D, R, r_exp = _parse_longest_route_large()
    r = find_longest_route_value_MT(n, m, D, R)
    assert r == r_exp


def test_longest_route_value_with_diagonal():
    D = np.array([[1, 1, 1], [1, 1, 1]])
    R = np.array([[1, 1], [1, 1], [1, 1]])
    DG = np.array([[1, 1], [1, 1]])
    assert find_longest_route_value_MTD(2, 2, D, R, DG) == 4
    D = np.array([[1, 1, 1], [1, 1, 1]])
    R = np.array([[4, 1], [1, 1], [1, 1]])
    DG = np.array([[1, 1], [1, 10]])
    assert find_longest_route_value_MTD(2, 2, D, R, DG) == 15


def test_construct():
    assert construct_route("AACC", "ACAC") == [
        ["s", "s", "s", "s", "s"],
        ["s", "dg", "r", "r", "r"],
        ["s", "d", "d", "dg", "r"],
        ["s", "d", "dg", "d", "dg"],
        ["s", "d", "d", "d", "d"],
    ]


def test_track_route():
    v = "AACC"
    w = "ACAC"
    assert track_longest_route(v, w) == "AAC"
    v = "AACCTTGG"
    w = "ACACTGTGA"
    assert track_longest_route(v, w) == "AACTTG"


def _parse_longest_common(file):
    with open(file) as fd:
        v = fd.readline().strip()
        w = fd.readline().strip()
    return v, w


TEST_DIR = Path("week1/data/LongestCommonSubsequence/inputs")


@pytest.mark.parametrize("inp", sorted(TEST_DIR.glob("*")))
def test_common_route(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    v, w = _parse_longest_common(inp)
    with open(out) as fd:
        cm_exp = fd.readline().strip()
    cm = track_longest_route(v, w)
    assert cm == cm_exp


def _parse_common_large():
    with open("week1/data/longest_common_subsequence.txt") as fd:
        fd.readline()
        v = fd.readline().strip()
        w = fd.readline().strip()
        fd.readline()
        cm = fd.readline().strip()
    return v, w, cm


def test_longest_common_large():
    v, w, cm_exp = _parse_common_large()
    cm = track_longest_route(v, w)
    assert cm == cm_exp


def test_convert_to_backtrack():
    graph = ((0, 1, 7), (1, 4, 1), (3, 4, 3))
    assert convert_to_backtrack(graph) == {1: [(0, 7)], 4: [(1, 1), (3, 3)]}


def test_find_backtrack():
    graph = ((0, 1, 7), (0, 2, 4), (2, 3, 2), (1, 4, 1), (3, 4, 3))
    assert find_highest_score_track(graph, 0, 4) == (9, [0, 2, 3, 4])


def _parse_backtrack_inp(file):
    with open(file) as fd:
        start = int(fd.readline().strip())
        end = int(fd.readline().strip())
        graph = []
        while True:
            line = fd.readline().strip()
            if not line:
                break
            line = re.split("->|:", line)
            graph.append([int(i.strip()) for i in line])
    return start, end, graph


def _parse_backtrack_out(file):
    with open(file) as fd:
        score = int(fd.readline().strip())
        L = fd.readline().strip().split("->")
        track = [int(i.strip()) for i in L]
    return score, track


def test_parse_backtrack(tmp_path):
    fn = tmp_path / "a.txt"
    fn.write_text(
        textwrap.dedent(
            """
    0
    10
    0->1:2
    3->4:5
    """
        ).strip()
    )
    st, end, gr = _parse_backtrack_inp(fn)
    assert st == 0
    assert end == 10
    assert gr == [[0, 1, 2], [3, 4, 5]]


TEST_DIR = Path("week1/data/LongestPathDAG/inputs/")


@pytest.mark.parametrize("inp", sorted(TEST_DIR.glob("*")))
def test_backtrack(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    start, end, graph = _parse_backtrack_inp(inp)
    score_exp, track_exp = _parse_backtrack_out(out)
    score, track = find_highest_score_track(graph, start, end)
    assert score == score_exp
    assert track == track_exp


def test_max_with_tie():
    assert max_with_tie(((1, 3), (2, 3)), 1) == [(1, 3), (2, 3)]


def _parse_longest_route_DAG_large(file):
    with open(file) as fd:
        fd.readline()
        start = int(fd.readline().strip())
        end = int(fd.readline().strip())
        graph = []
        while True:
            L = fd.readline().strip()
            if L == "":
                break
            L = re.split("->|:", L)
            graph.append([int(i.strip()) for i in L])
        fd.readline()
        score = int(fd.readline().strip())
        results = fd.readline().strip().split("->")
        path = [int(i.strip()) for i in results]
    return start, end, graph, score, path


def test_find_longest_route_DAG_large():
    start, end, graph, score, path = _parse_longest_route_DAG_large(
        "week1/data/longest_path_in_DAG.txt"
    )
    sc, pa = find_highest_score_track(graph, start, end)
    assert sc == score
    assert pa == path


def test_parse_longest_route_large(tmp_path):
    fn = tmp_path / "a.txt"
    fn.write_text(
        textwrap.dedent(
            """
    Input:
    0
    10
    0->1:2
    3->4:5

    Output:
    10
    0->1->2
    """
        ).strip()
    )
    st, end, gr, sc, path = _parse_longest_route_DAG_large(fn)
    assert st == 0
    assert end == 10
    assert gr == [[0, 1, 2], [3, 4, 5]]
    assert sc == 10
    assert path == [0, 1, 2]


def test_select_route_with_se():
    route = {0: (3, [0, 1, 2]), 1: (4, [1, 5])}
    assert select_route_with_startend(route, 0, 2) == (3, [0, 1, 2])
    assert select_route_with_startend(route, 1, 5) == (4, [1, 5])


def test_max_with_start_end():
    ls = [(0, 3), (1, 4), (2, 6)]
    assert max_with_start_end(ls, 0, 5) == [(0, 3)]
    assert max_with_start_end(ls, 1, 10) == [(1, 4)]
