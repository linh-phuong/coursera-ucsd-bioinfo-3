import numpy as np
from collections import defaultdict


def nonnegative_min(a_list):
    """find min non negative value of a list

    Args:
        a_list (list): a list

    Returns:
        float: min nonnegative value
    """
    ls = [i for i in a_list if i >= 0]
    return min(ls)


def min_nb_elements(S, elements):
    """find the number of minimum elements that add up to a sum

    Args:
        S (float): a number
        elements (list): a list of number

    Returns:
        int: number of minimum elements whose sum equals S
    """
    min_dict = {k: 1 for k in elements}
    elements = sorted(elements)
    for i in range(S + 1):
        ime = [i - e for e in elements]
        h = [min_dict.get(j, j) for j in ime]
        if i > elements[-1]:
            di = nonnegative_min(ime)
            min_dict.update({k: v for k, v in min_dict.items() if k >= di})
        m = nonnegative_min(h) + 1 if i > 0 else 0
        min_dict.update({i: m})
    return min_dict[S]


def find_min_elements(S, elements):
    """find the minimum elements that add up to a sum

    Args:
        S (float): a number
        elements (list): a list of numbers

    Returns:
        list: the elements that add up to the sum
    """
    min_dict = {k: [k] for k in elements}
    min_dict.update({0: []})
    max_e = max(elements)
    for i in range(1, S + 1):
        best, best_len = None, S + 1
        for e in elements:
            if e <= i:
                length = len(min_dict[i - e]) + 1
                if length < best_len:
                    best_len = length
                    best = min_dict[i - e] + [e]
        min_dict[i] = best
        if i > max_e:
            del min_dict[i - max_e - 1]
    return min_dict[S]


def find_longest_route_value_MT(c, r, D, R):
    """find the score of the longest route in a Manhattan graph (a rectangle graph in which each edge is a route)

    Args:
        c (int): number of vertical edges
        r (int): number of horizontal edges
        D (np.array): scores of vertical routes
        R (np.array): scores of horizontal routes

    Returns:
        int: score of the longest route from point (0, 0) to point (n, m) of the rectangle graph
    """
    n = c + 1
    m = r + 1
    route_score = np.zeros((n, m))
    for i in range(1, n):
        route_score[i, 0] = route_score[i - 1, 0] + D[i - 1, 0]
    for j in range(1, m):
        route_score[0, j] = route_score[0, j - 1] + R[0, j - 1]
    for i in range(1, n):
        for j in range(1, m):
            route_score[i, j] = max(
                (route_score[i - 1, j] + D[i - 1, j], route_score[i, j - 1] + R[i, j - 1])
            )
    return route_score[n - 1, m - 1]


def find_longest_route_value_MTD(c, r, D, R, DG):
    """find the score of the longest route in a Manhattan graph with diagonal route

    Args:
        c (int): number of vertical edges
        r (int): number of horizontal edges
        D (np.array): scores of vertical routes
        R (np.array): scores of horizontal routes
        DG (np.array): scores of diagonal routes

    Returns:
        int: score of the longest route from point (0, 0) to point (n, m) of the rectangle graph
    """
    n = c + 1
    m = r + 1
    route_score = np.zeros((n, m))
    for i in range(1, n):
        route_score[i, 0] = route_score[i - 1, 0] + D[i - 1, 0]
    for j in range(1, m):
        route_score[0, j] = route_score[0, j - 1] + R[0, j - 1]
    for i in range(1, n):
        for j in range(1, m):
            route_score[i, j] = max(
                (
                    route_score[i - 1, j] + D[i - 1, j],
                    route_score[i, j - 1] + R[i, j - 1],
                    route_score[i - 1, j - 1] + DG[i - 1, j - 1],
                )
            )
    return route_score[n - 1, m - 1]


def construct_route(v, w):
    """build directions to walk the longest route of a rectangle graph
    each column of the graph corresponds to each element of the string w
    each row of the graph corresponds to each element of the string v

    Args:
        v (str): a sequence of aa, or dna, etc
        w (str): a sequence of aa, or dna, etc

    Returns:
        list: directions from one node to the next node in the rectangle graph
        "s":stop, "d": go down, "r": go right, "dg": go diagonal
    """
    scan_v = len(v) + 1
    scan_w = len(w) + 1
    align = np.zeros((scan_v, scan_w))
    backtrack = [["s"] * scan_w]
    for i in range(1, scan_v):
        tr = ["s"]
        for j in range(1, scan_w):
            match = 1 if v[i - 1] == w[j - 1] else 0
            align[i, j] = max(align[i - 1, j], align[i, j - 1], align[i - 1, j - 1] + match)
            if align[i, j] == align[i - 1, j]:
                tr.append("d")
            elif align[i, j] == align[i, j - 1]:
                tr.append("r")
            elif align[i, j] == align[i - 1, j - 1] + match:
                tr.append("dg")
        backtrack.append(tr)
    return backtrack


def _track_longest_route(tracks, v, i, j):
    """find the longest route in a rectangle graph given the tracking information

    Args:
        tracks (list): directions from one node to the next node in the rectangle graph
        "s":stop, "d": go down, "r": go right, "dg": go diagonal
        v (int): number of vertical edges
        i, j (int, int): position (i, j) in the graph

    Returns:
        str: direction to the longest route from (0, 0) to (i, j)
    """
    if i == 0 or j == 0:
        return ""
    if tracks[i][j] == "d":
        return _track_longest_route(tracks, v, i - 1, j)
    elif tracks[i][j] == "r":
        return _track_longest_route(tracks, v, i, j - 1)
    else:
        return _track_longest_route(tracks, v, i - 1, j - 1) + v[i - 1]


def track_longest_route(v, w):
    """find the longest common string between v and w

    Args:
        v, w (str, str): a sequence of aa, dna, etc

    Returns:
        str: longest common string
    """
    tracks = construct_route(v, w)
    i = len(v)
    j = len(w)
    return _track_longest_route(tracks, v, i, j)


# graph: { 0: [(1, 3)] } => arc from 0 to 1 has weight 3
def find_highest_score_track(graph, start, end):
    """find the track from a start node to an end node with the highest score given a graph showing the weight when going from one node to another

    Args:
        graph (dict): showing edges from one node to other nodes with weight
        start (str or int): start node
        end (str or int): end node

    Raises:
        Exception: when the end node is not connected to the start node

    Returns:
        list: the path with the highest score
    """

    # backtrack: { 1: [(0, 3)] } => arc from 0 to 1 has weight 3
    backtrack = convert_to_backtrack(graph)

    allnodes = set(list(backtrack.keys()) + [n[0] for v in backtrack.values() for n in v])
    allnodes = sorted(allnodes)
    trimed_nodes = [i for i in allnodes if end >= i > start]

    # { 5: (7, [0, 3, 5])} => route from 0 to 5 with the highest score is 0-3-5 with score 7
    routes_dict = {start: (0, [start])}

    for i in trimed_nodes:
        bestroute, bestscore = None, -1
        for n, sc in backtrack[i]:  # edge from n to i has score sc
            if n in routes_dict:
                prev_score, prev_route = routes_dict[n]  # score and route from start to n
                score = sc + prev_score
                if score > bestscore:
                    bestroute, bestscore = prev_route + [i], score
        if bestroute is not None:
            routes_dict[i] = (bestscore, bestroute)
    if end not in routes_dict:
        raise Exception("The end node is not connected to the start node")
    return routes_dict[end]


def convert_to_backtrack(graph):
    """convert a graph in the list form, and each element is a list showing [from node, to node, weight] into
     a graph in the dictionary form where key = to node, value = from node and weight

    Args:
        graph (list): each element is a list [from node, to node, weight]

    Returns:
        defaultdict: key: to node, value = (weight, from node), defaut value is an empty list
    """
    backtrack = defaultdict(lambda: [])
    for items in graph:
        fr, to, weight = items
        backtrack[to].append((fr, weight))
    return backtrack
