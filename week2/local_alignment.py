def build_route_with_freeride(v, w, rpl):
    """build directions to walk the longest route of a rectangle graph
    with penalties for mismatch, insert and delete,
    and free ride from v[0], w[0] to any node v[j], w[i],
    a free ride means that accumulated punishment of mismatch and
    insert/delete is not taken into account
    each column of the graph corresponds to each element of the string w
    each row of the graph corresponds to each element of the string v

    Args:
        v (str): a sequence of aa, or dna, etc
        w (str): a sequence of aa, or dna, etc
        rpl (default_dict): matrix for match, mismatch,
        default value is the insert/deletion punishment

    Returns:
        list: directions from one node to the next node in the rectangle graph
        "d": go down, "r": go right, "dg": go diagonal
    """
    scan_v = len(v) + 1
    scan_w = len(w) + 1
    align = [[0 for j in range(scan_v)] for i in range(scan_w)]
    # if there is deletion or insertion, take the default value
    indel = rpl[()]
    backtrack = [["d"] * scan_v]
    maxscore = 0
    max_position = (0, 0)
    for i in range(1, scan_w):
        tr = ["r"]
        for j in range(1, scan_v):

            # if w[i-1] match or mismath with v[j-1]
            match_mismatch = rpl[(w[i - 1], v[j - 1])]

            # value at position [i, j] is the maximum from previous position
            # value = 0 means there is freeride from node [0, 0] to node [i, j]
            align[i][j] = max(
                0,
                align[i - 1][j] + indel,
                align[i][j - 1] + indel,
                align[i - 1][j - 1] + match_mismatch,
            )
            if align[i][j] > maxscore:
                maxscore = align[i][j]
                max_position = (i, j)
            if align[i][j] == align[i - 1][j] + indel:
                tr.append("d")
            elif align[i][j] == align[i][j - 1] + indel:
                tr.append("r")
            elif align[i][j] == align[i - 1][j - 1] + match_mismatch:
                tr.append("dg")
            else:
                tr.append("nc")
        backtrack.append(tr)
    return backtrack, maxscore, max_position


def match_local(v, w, i, j, backtrack):
    """Generate matched substring between two strings
    given a backtrack matrix and starting position [i, j]
    in the matrix

    Args:
        v (str): the first string
        w (str): the second string
        i (int): position on w
        j (int): position on v
        backtrack (list): back track matrix

    Returns:
        list: matched substrings between two strings
    """
    if i == 0 or j == 0 or backtrack[i][j] == "nc":
        return ["", ""]
    if backtrack[i][j] == "d":
        vv, ww = match_local(v, w, i - 1, j, backtrack)
        return [vv + "-", ww + w[i - 1]]
    elif backtrack[i][j] == "r":
        vv, ww = match_local(v, w, i, j - 1, backtrack)
        return [vv + v[j - 1], ww + "-"]
    else:
        vv, ww = match_local(v, w, i - 1, j - 1, backtrack)
        return [vv + v[j - 1], ww + w[i - 1]]


def local_alignment(v, w, replacement_score):
    """Create alignment strings between two strings

    Args:
        v (str): the first string
        w (str): the second string
        replacement_score (dict): score for match,
        punishment for mismatch and insert/deletion

    Returns:
        tuple: matched score, two matched substrings
    """
    backtrack, maxscore, max_position = build_route_with_freeride(v, w, replacement_score)
    wi, vj = max_position
    return maxscore, match_local(v, w, wi, vj, backtrack)