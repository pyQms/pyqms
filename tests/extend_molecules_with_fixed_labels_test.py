#!/usr/bin/env python
# encoding: utf-8


import pyqms


TESTS = [
    {
        "input": ["KLEINERTEST"],
        "params": {},
        "fixed_labels": {"R": ["Label:13C(6)15N(4)", ""]},
        "output": ["KLEINER0TEST", "KLEINER1TEST"],
    },
    {
        "input": ["KLEINERTEST"],
        "params": {},
        "fixed_labels": {
            "K": ["Label:13C(6)15N(2)", "Label:13C(6)15N(4)", ""],
            "R": ["Label:13C(6)15N(4)", ""],
        },
        "output": [
            "K0LEINER0TEST",
            "K0LEINER1TEST",
            "K1LEINER0TEST",
            "K1LEINER1TEST",
            "K2LEINER0TEST",
            "K2LEINER1TEST",
        ],
    },
    {
        "input": ["KLEINERTEST"],
        "params": {"SILAC_AAS_LOCKED_IN_EXPERIMENT": ["K", "R"]},
        "fixed_labels": {
            "K": ["Label:13C(6)15N(2)", ""],
            "R": ["Label:13C(6)15N(4)", ""],
        },
        "output": ["K0LEINER0TEST", "K1LEINER1TEST"],
    },
]


def generic_test():
    for test_dict in TESTS:
        yield generic_check_fucntion, test_dict


def generic_check_fucntion(test_dict):
    lib = pyqms.IsotopologueLibrary(
        molecules=test_dict["input"],
        charges=[2],
        params=test_dict["params"],
        fixed_labels=test_dict["fixed_labels"],
    )
    function_output = lib._extend_molecules_with_fixed_labels(test_dict["input"])
    print(function_output)
    # assert False
    assert sorted(list(function_output)) == sorted(test_dict["output"])


if __name__ == "__main__":
    for test_dict in TESTS:
        generic_check_fucntion(test_dict)
