DOWN = "d"
RIGHT = "r"
DIAGONAL = "dg"


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
    assert mu > 0, "input of punishment is positive"
    assert s > 0, "input of punishment is positive"
    # Initial score at each point
    align = [[0 for j in range(scan_v)] for i in range(scan_w)]

    # Initial direction to each point [i, j]
    backtrack = [[None for j in range(scan_v)] for i in range(scan_w)]

    for i in range(1, scan_w):
        align[i][0] = align[i - 1][0] - s
        backtrack[i][0] = DOWN
    for j in range(1, scan_v):
        align[0][j] = align[0][j - 1] - s
        backtrack[0][j] = RIGHT
    for i in range(1, scan_w):
        for j in range(1, scan_v):
            best = 0

            # calculate the weight when going
            # down, right, and diagonal to reach [i, j]
            wt_r = align[i][j - 1] - s
            wt_d = align[i - 1][j] - s
            wt_dg = align[i - 1][j - 1] + m if w[i - 1] == v[j - 1] else align[i - 1][j - 1] - mu

            # value at position [i, j] is the maximum from previous position
            best = max(wt_d, wt_r, wt_dg)
            align[i][j] = best

            if best == wt_d:
                backtrack[i][j] = DOWN
            elif best == wt_r:
                backtrack[i][j] = RIGHT
            elif best == wt_dg:
                backtrack[i][j] = DIAGONAL
            else:
                raise Exception("No connection with previous position")
    return backtrack, align[-1][-1]


# weight = {("A", "A"): 3, ("A", "B"):-2}
# A -> A: matched, plus 3 point, A -> B: mismatch, minus 2 point
# weight[()] = default value for insert/deletion of an element
def build_route_replace_matrix(v, w, weight):
    """build directions to walk the longest route of a rectangle graph
    with penalties for mismatch, insert and delete
    each column of the graph corresponds to each element of the string w
    each row of the graph corresponds to each element of the string v

    Args:
        v (str): a sequence of aa, or dna, etc
        w (str): a sequence of aa, or dna, etc
        weight (default_dict): score for match, mismatch, and
        insert/deletion

    Returns:
        list: directions from one node to the next node in the rectangle graph
        "d": go down, "r": go right, "dg": go diagonal
    """
    scan_v = len(v) + 1
    scan_w = len(w) + 1

    # Initial score at each point
    align = [[0 for j in range(scan_v)] for i in range(scan_w)]

    # Initial direction to each point [i, j]
    backtrack = [[None for j in range(scan_v)] for i in range(scan_w)]

    # default value of the weight matrix
    # is punishment for insertion/deletion
    indel = weight[()]

    for i in range(1, scan_w):
        align[i][0] = align[i - 1][0] + indel
        backtrack[i][0] = DOWN
    for j in range(1, scan_v):
        align[0][j] = align[0][j - 1] + indel
        backtrack[0][j] = RIGHT
    for i in range(1, scan_w):
        for j in range(1, scan_v):
            best = 0
            # if w[i-1] match or mismath with v[j-1]
            match_mismatch = weight[(w[i - 1], v[j - 1])]

            # calculate the weight when going
            # down, right, and diagonal to reach [i, j]
            wt_r = align[i][j - 1] + indel
            wt_d = align[i - 1][j] + indel
            wt_dg = align[i - 1][j - 1] + match_mismatch

            # value at position [i, j] is the maximum from previous position
            best = max(wt_d, wt_r, wt_dg)
            align[i][j] = best

            if best == wt_d:
                backtrack[i][j] = DOWN
            elif best == wt_r:
                backtrack[i][j] = RIGHT
            elif best == wt_dg:
                backtrack[i][j] = DIAGONAL
            else:
                raise Exception("No connection with previous position")
    return backtrack, align[-1][-1]


def match_strings(v, w, backtrack):
    """Globally align two strings given a backtrack matrix

    Args:
        v (str): the first string
        w (str): the second string
        backtrack (list): backtrack matrix

    Returns:
        tuple: len = 2, global alignment of two strings
    """
    i = len(w)
    j = len(v)
    v_align, w_align = "", ""
    if backtrack[i][j] == None:
        if i == 0 and j == 0:
            return None
        else:
            raise Exception("The current position does not connect to any previous position")
    while i != 0 or j != 0:
        if backtrack[i][j] == "d":
            v_align += "-"
            w_align += w[i - 1]
            i -= 1
        elif backtrack[i][j] == "r":
            v_align += v[j - 1]
            w_align += "-"
            j -= 1
        else:
            v_align += v[j - 1]
            w_align += w[i - 1]
            i -= 1
            j -= 1
    return v_align[::-1], w_align[::-1]


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
    return score, match_strings(v, w, backtrack)


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
    return score, match_strings(v, w, backtrack)