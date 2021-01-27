from week1.common_element import find_longest_common, find_min_elements, min_nb_elements


def test_find_longest_common():
    assert find_longest_common("ACTGCA", "CATCGT") == "ACT"
    assert find_longest_common("ACTGCA", "CATCGC") == "ACGC"
    assert find_longest_common("AACTG", "CAAG") == "AAG"
    assert find_longest_common("ATGTTATA", "ATCGTCC") == "ATGT"


def test_min_nb_elements():
    c = [1, 4, 5]
    assert min_nb_elements(0, c) == 0
    assert min_nb_elements(1, c) == 1
    assert min_nb_elements(2, c) == 2
    assert min_nb_elements(3, c) == 3
    assert min_nb_elements(4, c) == 1
    assert min_nb_elements(22, c) == 5
    c = [50, 25, 20, 10, 5, 1]
    assert min_nb_elements(40, c) == 2


def test_find_min_elements():
    c = [1, 4, 5]
    assert find_min_elements(0, c) == []
    assert find_min_elements(1, c) == [1]
    assert find_min_elements(2, c) == [1, 1]
    assert find_min_elements(4, c) == [4]
    assert find_min_elements(22, c) == [5, 5, 4, 4, 4]
