from week2.local_alignment import match_local
from week2.global_alignment import global_aligment, global_alignment_repmat


def hamming_distance(v, w):
    assert len(v) == len(w), "lengths of two strings are not equal"
    d = 0
    for i, j in zip(v, w):
        if i != j:
            d += 1
    return d


def calculate_edit_distance(v, w, penalties):
    _, align = global_aligment(v, w, penalties)
    hd = hamming_distance(*align)
    return hd


def calculate_edit_distance_replacemat(v, w, replace_matrix):
    _, align = global_alignment_repmat(v, w, replace_matrix)
    hd = hamming_distance(*align)
    return hd


DOWN = "d"
RIGHT = "r"
DIAGONAL = "dg"
FREERIDE = "fr"


def build_route_with_freeride_to_first_match(v, w, weight):
    scan_v = len(v) + 1
    scan_w = len(w) + 1

    # default value of the weight matrix
    # is punishment for insertion/deletion
    indel = weight[()]
    # Initial score at each point
    align = [[0 for j in range(scan_v)] for i in range(scan_w)]

    # Initial direction to each point [i, j]
    backtrack = [[None for j in range(scan_v)] for i in range(scan_w)]

    # Provide free ride from starting point of v to first matching point
    for j in range(scan_v):
        backtrack[0][j] = FREERIDE
    for i in range(1, scan_w):
        backtrack[i][0] = DOWN
        align[i][0] = align[i - 1][0] + indel

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

            # No free ride if it is in the middle of w
            # but free ride if the scan on string w is finished
            if i < scan_w:
                best = max(wt_d, wt_r, wt_dg)
            else:
                best = max(0, wt_d, wt_r, wt_dg)

            align[i][j] = best

            if best == wt_d:
                backtrack[i][j] = DOWN
            elif best == wt_r:
                backtrack[i][j] = RIGHT
            else:
                backtrack[i][j] = DIAGONAL
    return backtrack, align


def find_fitting_positions(align_score):
    max_score = max(align_score[-1])
    last_row_id = len(align_score) - 1
    positions = []
    for j, sc in enumerate(align_score[-1]):
        # the position has to be at least the last point of w
        if sc == max_score and j > last_row_id - 1:
            positions.append((last_row_id, j))
    return max_score, positions


def find_fitting_alignment(v, w, weight):
    """Find alignment that fit a string (w) to another string (v)

    Args:
        v (str): the first string
        w (str): the second string
        weight (list): matrix with score when there is a match,
        punishment when there is a mismatch and insert/deletion

    Returns:
        tuple: score of alignment, list of possible aligments
    """
    bt, al_score = build_route_with_freeride_to_first_match(v, w, weight)
    max_score, positions = find_fitting_positions(al_score)
    alignments = []
    for i, j in positions:
        _, a = match_local(v, w, i, j, bt)
        alignments.append(match_local(v, w, i, j, bt))
    return max_score, alignments


def find_overlap_positions(align_score):
    last_col_id = len(align_score[0]) - 1
    len_rows = len(align_score)
    best_score = -1
    positions = []
    for i in range(1, len_rows):
        if align_score[i][last_col_id] > best_score:
            best_score = align_score[i][last_col_id]
            if positions:
                positions = [(i, last_col_id)]
            else:
                positions.append((i, last_col_id))
    return best_score, positions


def find_overlap_alignment(v, w, weight):
    bt, al_score = build_route_with_freeride_to_first_match(v, w, weight)
    max_score, positions = find_overlap_positions(al_score)
    alignments = []
    for i, j in positions:
        _, a = match_local(v, w, i, j, bt)
        alignments.append(match_local(v, w, i, j, bt))
    return max_score, alignments