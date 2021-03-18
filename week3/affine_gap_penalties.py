from week2.global_alignment import match_strings
import numpy as np

INSERT = "d"
DELETE = "r"
DIAGONAL = "dg"


def build_graph_with_affine_gap_pen(v, w, weights):
    """build backtrack with affine gap penalties

    Args:
        v (str): the first string
        w (str): the second string
        weights (tuple): weights for each movement within and between graphs

    Returns:
        tuple: score of the best route and the backtrack matrix
    """

    # replacement is a dictionary for match/mismatch point, s equals to an insertion or deletion,
    # and e equals to insertion or deletion with affine gap
    replacement, s, e = weights
    col = len(w) + 1
    row = len(v) + 1

    # Default values of the graph are negative infinity
    diag_graph = [[-np.inf for j in range(col)] for i in range(row)]
    insert_graph = [[-np.inf for j in range(col)] for i in range(row)]
    delete_graph = [[-np.inf for j in range(col)] for i in range(row)]
    direction_graph = [[None for j in range(col)] for i in range(row)]

    # Build values for the fist row and first column

    diag_graph[0][0] = 0
    for j in range(1, col):
        delete_graph[0][j] = max(delete_graph[0][j - 1] - e, diag_graph[0][j - 1] - s)
        diag_graph[0][j] = max(insert_graph[0][j], delete_graph[0][j])
        direction_graph[0][j] = DELETE
    for i in range(1, row):
        insert_graph[i][0] = max(insert_graph[i - 1][0] - e, diag_graph[i - 1][0] - s)
        diag_graph[i][0] = max(insert_graph[i][0], delete_graph[i][0])
        direction_graph[i][0] = INSERT

    # Build values for the point [i][j] based on previous point values

    for j in range(1, col):
        for i in range(1, row):
            sc = replacement[(v[i - 1], w[j - 1])]
            insert_graph[i][j] = max(insert_graph[i - 1][j] - e, diag_graph[i - 1][j] - s)
            delete_graph[i][j] = max(delete_graph[i][j - 1] - e, diag_graph[i][j - 1] - s)
            diag_graph[i][j] = max(
                delete_graph[i][j], insert_graph[i][j], diag_graph[i - 1][j - 1] + sc
            )
            best_score = max(insert_graph[i][j], delete_graph[i][j], diag_graph[i][j])
            if best_score == insert_graph[i][j]:
                direction_graph[i][j] = INSERT
            elif best_score == delete_graph[i][j]:
                direction_graph[i][j] = DELETE
            else:
                direction_graph[i][j] = DIAGONAL
    accum_score = max(
        diag_graph[row - 1][col - 1], insert_graph[row - 1][col - 1], delete_graph[row - 1][col - 1]
    )
    return direction_graph, accum_score


def align_with_affine_gap_pen(v, w, weights):
    """align two strings with affine gap penalties

    Args:
        v (str): string one
        w (str): string two
        weights (tuple): weights for movement within and between graphs

    Returns:
        tuple: score of the best alignment, alignment
    """
    backtrack, score = build_graph_with_affine_gap_pen(v, w, weights)
    # in the match_strings function w is the first string and v is the second string
    alignment = match_strings(w, v, backtrack)
    return score, alignment[::-1]