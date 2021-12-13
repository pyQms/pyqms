#!/usr/bin/env python
# encoding: utf-8
"""

Test _extend_kb_with_fixed_labels from core

"""
import pyqms
import sys
import unittest

TESTS = [
    {
        "in": {
            "params": {
                "molecules": ["UU"],
                "charges": [2],
                "params": {
                    "AMINO_ACIDS": {"U": "C(3)H(7)N(1)O(2)Se(1)"},
                },
            }
        },
        "out": {
            "aa": {"U": {"C": 3, "H": 7, "N": 1, "O": 2, "Se": 1}},
            "formula": "C(6)H(16)N(2)O(5)Se(2)"
        },
    },
]
# 2 isotope element (N,nitrogen)
CRASH_TESTS = {
    "in": {
        "params": {
            "molecules": ["KLEINERTEST"],
            "charges": [2],
            "fixed_labels": {
                # non existing aa
                "X": ["C(-6) 13C(6) N(-4) 15N(4)"]
            },
        }
    },
    "out": {
        # 'formated_molecule' : ['KLEINER0TEST'],
    },
}


def extend_kb_with_amino_acids_test():
    for test_id, test_dict in enumerate(TESTS):
        yield _extend_kb_with_amino_acids, test_id, test_dict


def _extend_kb_with_amino_acids(test_id, test_dict):

    lib_1 = pyqms.IsotopologueLibrary(**test_dict["in"]["params"])
    for aa, composition in test_dict["out"]["aa"].items():
        assert aa in lib_1.aa_compositions
        print(f"lib1 composition: {str(lib_1.aa_compositions[aa])}")
        print(f"expected composition: {composition}")
        assert composition == lib_1.aa_compositions[aa]

    formula_1 = list(lib_1.keys())[0]
    assert test_dict["out"]["formula"] == formula_1


class TestResults(unittest.TestCase):
    def setUp(self):
        pass

    def crash_test(self):
        """
        Check if a key error is raised when using a non existent amino acid
        """
        with self.assertRaises(SystemExit) as system_exit_check:
            pyqms.IsotopologueLibrary(**CRASH_TESTS["in"]["params"])
        self.assertEqual(system_exit_check.exception.code, 1)


if __name__ == "__main__":
    pass
