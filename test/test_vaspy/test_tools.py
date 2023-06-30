#

from vaspy import tools


def test_atom_selection_to_list():
    assert tools.atom_selection_to_list("1-5,8,8,9-15,10", False) == [
        "1",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "2",
        "3",
        "4",
        "5",
        "8",
        "9",
    ]
    assert tools.atom_selection_to_list("1-5,8,8,9-15,10") == [
        1,
        2,
        3,
        4,
        5,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
    ]
