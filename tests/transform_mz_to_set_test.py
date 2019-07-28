#!/usr/bin/env python
# encoding: utf-8
"""

Testfunction to test _transform_mz_to_set from core

"""

import pyqms


TESTS = {
    "mo1": {"i": 1000.001, "o": 21, "kwargs": None, "params": {"REL_MZ_RANGE": 10e-6}},
    "mo2": {
        "i": 1000.0009,
        "o": 101,
        "kwargs": None,
        "params": {"REL_MZ_RANGE": 5e-6, "INTERNAL_PRECISION": 10000},
    },
}


def transformation_mz_test():
    for test_id, test_dict in TESTS.items():
        yield transformation_mz, test_id, test_dict


def transformation_mz(test_id, test_dict):
    lib = pyqms.IsotopologueLibrary(molecules=["KLEINERTEST"], charges=[2])
    lib.params.update(test_dict["params"])
    tmp = lib._transform_mz_to_set(test_dict["i"])
    assert len(tmp) == test_dict["o"]


if __name__ == "__main__":
    transformation_mz_test()
