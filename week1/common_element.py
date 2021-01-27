def find_longest_common(a, b):
    """find the longest common string between a and b (duplicates and position included)

    Args:
        a (str): string a
        b (str): string b

    Returns:
        str: common string between a and b
    """
    c = ""
    i = 0
    for e in a:
        new_b = b[i:]
        j = 0
        if not new_b:
            break
        for f in new_b:
            i += 1
            if e == f:
                c += e
                break
            else:
                j += 1
            if i == len(b):
                i = i - j
    return c


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
