#!/usr/bin/env python
# encoding: utf-8
"""

Test _extend_kb_with_fixed_labels from core

"""
import pyqms
import sys
import unittest

TESTS = [
    # {
    #     'in' : {
    #         'params' : {
    #             'molecules' : ['KLEINERTEST'],
    #             'charges' : [2, ],
    #             'fixed_labels' : {
    #                 'R' : ['C(-6) 13C(6) N(-4) 15N(4)']
    #             },
    #         }
    #     },
    #     'out' : {
    #         'formated_molecule' : ['KLEINER0TEST'],
    #     }
    # },
    {
        "in": {
            "params": {
                "molecules": ["KLEINERTEST"],
                "charges": [2],
                "fixed_labels": {
                    "R": ["", "C(-6) 13C(6) N(-4) 15N(4)"],
                    "K": ["", "C(-6) 13C(6) N(-4) 15N(4)"],
                },
            }
        },
        "out": {
            "formated_molecule": sorted(
                ["K0LEINER0TEST", "K1LEINER0TEST", "K1LEINER1TEST", "K0LEINER1TEST"]
            )
        },
    },
    {
        "in": {
            "params": {
                "molecules": ["KLEINERTEST"],
                "charges": [2],
                "fixed_labels": {
                    "R": ["", "C(-6) 13C(6) N(-4) 15N(4)"],
                    "K": ["", "C(-6) 13C(6) N(-4) 15N(4)"],
                },
                "params": {"SILAC_AAS_LOCKED_IN_EXPERIMENT": ["K", "R"]},
            }
        },
        "out": {"formated_molecule": sorted(["K0LEINER0TEST", "K1LEINER1TEST"])},
    },
    {
        "in": {
            "params": {
                "molecules": ["KLEINERTEST"],
                "charges": [2],
                "fixed_labels": {
                    "R": ["", "C(-6) 13C(6) N(-4) 15N(4)"],
                    "K": ["", "C(-6) 13C(6) N(-4) 15N(4)"],
                    "I": ["FOCK"],
                },
                "params": {"SILAC_AAS_LOCKED_IN_EXPERIMENT": ["K", "R"]},
            }
        },
        "out": {"formated_molecule": sorted(["K0LEI0NER0TEST", "K1LEI0NER1TEST"])},
    },
    {
        "in": {
            "params": {
                "molecules": ["KLEINERTEST"],
                "charges": [2],
                "fixed_labels": {
                    "R": ["", "C(-6) 13C(6) N(-4) 15N(4)"],
                    "K": ["", "C(-6) 13C(6) N(-4) 15N(4)"],
                    "I": ["FOCK", ""],
                },
                "params": {"SILAC_AAS_LOCKED_IN_EXPERIMENT": ["K", "R"]},
            }
        },
        "out": {
            "formated_molecule": sorted(
                ["K0LEI0NER0TEST", "K1LEI0NER1TEST", "K0LEI1NER0TEST", "K1LEI1NER1TEST"]
            )
        },
    },
    {
        "in": {
            "params": {
                "molecules": ["KLEINERTEST"],
                "charges": [2],
                "fixed_labels": {
                    "R": ["", "C(-6) 13C(6) N(-4) 15N(4)"],
                    "K": ["", "C(-6) 13C(6) N(-4) 15N(4)"],
                    "I": ["FOCK", ""],
                    "L": ["Noo", "Way"],
                },
                "params": {"SILAC_AAS_LOCKED_IN_EXPERIMENT": ["K", "R"]},
            }
        },
        "out": {
            "formated_molecule": sorted(
                [
                    "K0L0EI0NER0TEST",
                    "K1L0EI0NER1TEST",
                    "K0L0EI1NER0TEST",
                    "K1L0EI1NER1TEST",
                    "K0L1EI0NER0TEST",
                    "K1L1EI0NER1TEST",
                    "K0L1EI1NER0TEST",
                    "K1L1EI1NER1TEST",
                ]
            )
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
                "U": ["C(-6) 13C(6) N(-4) 15N(4)"]
            },
        }
    },
    "out": {
        # 'formated_molecule' : ['KLEINER0TEST'],
    },
}


def extend_kb_with_fixed_labels_test():
    for test_id, test_dict in enumerate(TESTS):
        yield _extend_kb_with_fixed_labels, test_id, test_dict


def _extend_kb_with_fixed_labels(test_id, test_dict):

    lib_1 = pyqms.IsotopologueLibrary(**test_dict["in"]["params"])
    print(lib_1.lookup["molecule fixed label variations"])
    formula_1 = list(lib_1.keys())[0]
    # __oOo__
    lookup_key = test_dict["in"]["params"]["molecules"][0]
    for label_percentile in lib_1[formula_1]["env"].keys():
        assert (
            sorted(list(lib_1.lookup["molecule fixed label variations"][lookup_key]))
            == test_dict["out"]["formated_molecule"]
        )


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
