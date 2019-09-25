#!/usr/bin/env python
# encoding: utf-8
"""

Testfunction to test _transform_mz_to_set from core

"""

import pyqms


TESTS = {
    "test1": {
        "i": [(1000.0009, 1)],
        "o": {"tmz_set": set([10000009]), "tmz_lookup": {10000009: [(1000.0009, 1)]}},
        "kwargs": None,
        "params": {"REL_MZ_RANGE": 5e-6, "INTERNAL_PRECISION": 10000},
    }
}


def transformation_spectrum_test():
    for test_id, test_dict in TESTS.items():
        yield transformation_spectrum, test_id, test_dict


def transformation_spectrum(test_id, test_dict):
    lib = pyqms.IsotopologueLibrary(molecules=["KLEINERTEST"], charges=[2])
    lib.params.update(test_dict["params"])
    tmz_set, tmz_lookup = lib._transform_spectrum(test_dict["i"], mz_range=None)
    assert sorted(tmz_set) == sorted(test_dict["o"]["tmz_set"])
    for tmz in tmz_set:
        assert tmz in tmz_lookup.keys()
        assert test_dict["o"]["tmz_lookup"][tmz] == test_dict["i"]

    for tmz in tmz_lookup.keys():
        assert tmz_lookup[tmz] == test_dict["o"]["tmz_lookup"][tmz]


if __name__ == "__main__":
    transformation_spectrum_test()
