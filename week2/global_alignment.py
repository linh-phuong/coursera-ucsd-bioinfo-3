def build_route_with_pen(v, w, penalties):
    """build directions to walk the longest route of a rectangle graph
    with penalties for mismatch, insert and delete
    each column of the graph corresponds to each element of the string w
    each row of the graph corresponds to each element of the string v

    Args:
        v (str): a sequence of aa, or dna, etc
        w (str): a sequence of aa, or dna, etc
        penalties (tuple): score if match, mismatch and insert/delete

    Returns:
        list: directions from one node to the next node in the rectangle graph
        "d": go down, "r": go right, "dg": go diagonal
    """
    scan_v = len(v) + 1
    scan_w = len(w) + 1
    m, mu, s = penalties
    align = [[0 for j in range(scan_v)] for i in range(scan_w)]
    for i in range(1, scan_w):
        align[i][0] = align[i - 1][0] - s
    for j in range(1, scan_v):
        align[0][j] = align[0][j - 1] - s
    backtrack = [["d"] * scan_v]
    for i in range(1, scan_w):
        tr = ["r"]
        for j in range(1, scan_v):
            match_mismatch = m if v[j - 1] == w[i - 1] else -mu
            align[i][j] = max(
                align[i - 1][j] - s, align[i][j - 1] - s, align[i - 1][j - 1] + match_mismatch
            )
            if align[i][j] == align[i - 1][j] - s:
                tr.append("d")
            elif align[i][j] == align[i][j - 1] - s:
                tr.append("r")
            elif align[i][j] == align[i - 1][j - 1] + match_mismatch:
                tr.append("dg")
            else:
                raise Exception("No connection with previous position")
        backtrack.append(tr)
    return backtrack, align[-1][-1]


def build_route_replace_matrix(v, w, rpl):
    """build directions to walk the longest route of a rectangle graph
    with penalties for mismatch, insert and delete
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

    for i in range(1, scan_w):
        align[i][0] = align[i - 1][0] + indel
    for j in range(1, scan_v):
        align[0][j] = align[0][j - 1] + indel
    backtrack = [["d"] * scan_v]
    for i in range(1, scan_w):
        tr = ["r"]
        for j in range(1, scan_v):

            # if w[i-1] match or mismath with v[j-1]
            match_mismatch = rpl[(w[i - 1], v[j - 1])]

            # value at position [i, j] is the maximum from previous position
            align[i][j] = max(
                align[i - 1][j] + indel,
                align[i][j - 1] + indel,
                align[i - 1][j - 1] + match_mismatch,
            )
            if align[i][j] == align[i - 1][j] + indel:
                tr.append("d")
            elif align[i][j] == align[i][j - 1] + indel:
                tr.append("r")
            elif align[i][j] == align[i - 1][j - 1] + match_mismatch:
                tr.append("dg")
            else:
                raise Exception("No connection with previous position")
        backtrack.append(tr)
    return backtrack, align[-1][-1]


def match_strings(v, w, i, j, backtrack):
    """Globally align two strings given a backtrack matrix
    and position [i, j] in the matrix

    Args:
        v (str): the first string
        w (str): the second string
        i (int): position on w
        j (int): position on v
        backtrack (list): backtrack matrix

    Returns:
        list: len = 2, global alignment of two strings
    """
    if i == 0 and j == 0:
        return ["", ""]
    elif i == 0 and j > 0:
        vv, ww = match_strings(v, w, i, j - 1, backtrack)
        return [vv + v[j - 1], ww + "-"]
    elif i > 0 and j == 0:
        vv, ww = match_strings(v, w, i - 1, j, backtrack)
        return [vv + "-", ww + w[i - 1]]
    if backtrack[i][j] == "d":
        vv, ww = match_strings(v, w, i - 1, j, backtrack)
        return [vv + "-", ww + w[i - 1]]
    elif backtrack[i][j] == "r":
        vv, ww = match_strings(v, w, i, j - 1, backtrack)
        return [vv + v[j - 1], ww + "-"]
    else:
        vv, ww = match_strings(v, w, i - 1, j - 1, backtrack)
        return [vv + v[j - 1], ww + w[i - 1]]


def global_aligment(v, w, penalties):
    """Create global alignment given
    score for match and penalties for mismatch and insert/deletion

    Args:
        v (str): the first string
        w (str): the second string
        penalties (list): match, mismatch, insert/deletion

    Returns:
        tuple: matched score, alignments
    """
    backtrack, score = build_route_with_pen(v, w, penalties)
    lv, lw = len(v), len(w)
    return score, match_strings(v, w, lw, lv, backtrack)


def global_alignment_repmat(v, w, replacement_score):
    """Create global alignment given
    score for match and penalties for mismatch and insert/deletion

    Args:
        v (str): the first string
        w (str): the second string
        replacement_score (dict): score for match, punishment for mismatch and insert/deletion

    Returns:
        tuple: matched score, alignments
    """
    backtrack, score = build_route_replace_matrix(v, w, replacement_score)
    lv, lw = len(v), len(w)
    return score, match_strings(v, w, lw, lv, backtrack)