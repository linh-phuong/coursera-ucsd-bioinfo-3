def hamming_distance(v, w):
    d = 0
    for i, j in zip(v, w):
        if i != j:
            d += 1
    return d
