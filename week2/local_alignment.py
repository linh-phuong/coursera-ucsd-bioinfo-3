DOWN = "d"
RIGHT = "r"
DIAGONAL = "dg"
FREERIDE = "fr"


def build_route_with_freeride(v, w, weight):
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
        list of maximum score, position of the maximum score
        and matrix showing directions from one node to the next node in the rectangle graph
        "d": go down, "r": go right, "dg": go diagonal
    """
    scan_v = len(v) + 1
    scan_w = len(w) + 1
    # Initial score at each point
    align = [[0 for j in range(scan_v)] for i in range(scan_w)]

    # Initial direction to each point [i, j]
    backtrack = [[FREERIDE for j in range(scan_v)] for i in range(scan_w)]

    # default value of the weight matrix
    # is punishment for insertion/deletion
    indel = weight[()]
    maxscore = 0
    max_position = (0, 0)
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

            # value at position [i, j] is the maximum from previous positions
            # above, on the left, on the diagonal or from the souce [0, 0]
            best = max(wt_d, wt_r, wt_dg, 0)
            align[i][j] = best
            if best > maxscore:
                maxscore = best
                max_position = (i, j)

            if best == wt_d:
                backtrack[i][j] = DOWN
            elif best == wt_r:
                backtrack[i][j] = RIGHT
            elif best == wt_dg:
                backtrack[i][j] = DIAGONAL
            else:
                backtrack[i][j] = FREERIDE
    return backtrack, maxscore, max_position


def match_local(v, w, i, j, backtrack):
    """Align two strings locally given a starting position (i, j)
    and a backtrack matrix

    Args:
        v (str): the first string
        w (str): the second string
        i (int): row position in the backtrack matrix
        j (int): collum position in the backtrack matrix
        backtrack (list): give direction from [i, j] to a freeride point from the source [0, 0]

    Raises:
        Exception: when [i, j] does not connect to any position

    Returns:
        tuple: aligned substrings of v and w
    """
    v_align, w_align = "", ""
    if backtrack[i][j] is None:
        if i == 0 and j == 0:
            return None
        else:
            raise Exception("The current position does not connect to any previous position")
    while i != 0 or j != 0:
        if backtrack[i][j] == FREERIDE:
            break
        if backtrack[i][j] == DOWN:
            v_align += "-"
            w_align += w[i - 1]
            i -= 1
        elif backtrack[i][j] == RIGHT:
            v_align += v[j - 1]
            w_align += "-"
            j -= 1
        else:
            v_align += v[j - 1]
            w_align += w[i - 1]
            i -= 1
            j -= 1
    return v_align[::-1], w_align[::-1]


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
    backtrack, maxscore, max_pos = build_route_with_freeride(v, w, replacement_score)
    wi, vj = max_pos
    return maxscore, match_local(v, w, wi, vj, backtrack)
