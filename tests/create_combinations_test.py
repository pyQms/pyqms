#!/usr/bin/env python
# encoding: utf-8


import pyqms

lib = pyqms.IsotopologueLibrary(molecules=["PAINLESS"], charges=[2])

TESTS = [
    {
        "input": [("X", 2), ("Y", 2)],
        "output": [
            [("X", 0), ("Y", 0)],
            [("X", 0), ("Y", 1)],
            [("X", 1), ("Y", 0)],
            [("X", 1), ("Y", 1)],
        ],
    },
    {
        "input": [("X", 1), ("Y", 1), ("Z", 1)],
        "output": [[("X", 0), ("Y", 0), ("Z", 0)]],
    },
    {
        "input": [("X", 1), ("Y", 2), ("Z", 1)],
        "output": [[("X", 0), ("Y", 0), ("Z", 0)], [("X", 0), ("Y", 1), ("Z", 0)]],
    },
]


def generic_test():
    for test_dict in TESTS:
        yield generic_check_fucntion, test_dict


def generic_check_fucntion(test_dict):
    function_output = lib._create_combinations(
        test_dict["input"], all_combos=None, current_combo=None, pos=0
    )
    print(function_output)
    assert sorted(function_output) == test_dict["output"]


if __name__ == "__main__":
    for test_dict in TESTS.keys():
        generic_check_fucntion(test_dict)
