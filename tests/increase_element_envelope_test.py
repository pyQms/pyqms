#!/usr/bin/env python
# encoding: utf-8


import pyqms


TESTS = [
    {
        "input": ["KLEINERTEST"],
        "metabolic_labels": {"18O": [0, 0.99]},
        "output": [(("O", "0.000"),), (("O", "0.990"),)],
    }
]


def generic_test():
    for test_dict in TESTS:
        yield generic_check_fucntion, test_dict


def generic_check_fucntion(test_dict):
    lib = pyqms.IsotopologueLibrary(
        molecules=test_dict["input"],
        charges=[2],
        metabolic_labels=test_dict["metabolic_labels"],
    )
    for tuple_2_check in lib.labled_percentiles:
        assert tuple_2_check in test_dict["output"]
    assert len(lib.labled_percentiles) == len(test_dict["output"])


if __name__ == "__main__":
    for test_dict in TESTS:
        generic_check_fucntion(test_dict)
