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

def build_fitting_route(v, w, weight):
    
