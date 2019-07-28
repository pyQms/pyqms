#!/usr/bin/env python
# encoding: utf-8
"""
 Example of data passed::
    {
        'params': {'measurement_and_reporting': {'NAME': 'default'}},
        'fixed_labels': [
            {
                'modification': {
                    'unimodID': '4',
                    'specificity_sites': ['C'],
                    'mono_mass': 57.021464,
                    'element': {'O': 1, 'H': 3, 'N': 1, 'C': 2},
                    'name': 'Carbamidomethyl'
                },
                'AA': 'C'
            }
        ],
        'molecules': 'AA',
        'metabolic_labels': [{'modification': '0, 0.99', 'atom': '15N'}],
        'charges': [1, 2, 3, 4, 5],
        'file': '/BSA1.mzML'
    }

"""

from pyqms.adaptors import (
    _parse_evidence_and_format_fixed_labels as _parse_evidence_and_format_fixed_labels,
)
import os


TESTS = [
    {
        "input": {
            "params": {"measurement_and_reporting": {"NAME": "default"}},
            "fixed_labels": [
                {
                    "modification": {
                        "unimodID": "4",
                        "specificity_sites": ["C"],
                        "mono_mass": 57.021464,
                        "element": {"O": 1, "H": 3, "N": 1, "C": 2},
                        "name": "Carbamidomethyl",
                    },
                    "AA": "C",
                },
                {"modification": {"name": "None"}, "AA": "C"},
            ],
            "charges": [],
            "molecules": [
                "CLEINERTEST#Carbamidomethyl:1",
                "CLEINERTEST#Carbamidomethyl:1;Oxidation:2",
            ],
            "evidence_score_field": "PEP",
        },
        "output": {
            "fixed_labels": {"C": ["C(2)H(3)14N(1)O(1)", ""]},
            "molecules": ["CLEINERTEST", "CLEINERTEST#Oxidation:2"],
        },
    },
    {
        "input": {
            "params": {"measurement_and_reporting": {"NAME": "default"}},
            "charges": ["True", "False", "2", "3", "4", "5"],
            "molecules": [],
            "evidence_score_field": "PEP",
        },
        "output": {"charges": set([2, 3, 4, 5])},
    },
    {
        "input": {
            "params": {
                "measurement_and_reporting": {
                    "NAME": "default",
                    "MACHINE_OFFSET_IN_PPM": "5",
                    "MZ_SCORE_PERCENTILE": "0.6",
                }
            },
            "charges": [],
            "molecules": [],
            "evidence_score_field": "PEP",
        },
        "output": {
            "params": {
                # 'NAME'                  : 'default',
                "MACHINE_OFFSET_IN_PPM": 5.0,
                "MZ_SCORE_PERCENTILE": 0.6,
            }
        },
    },
    {
        # test to cover some execptions and else blocks...
        "input": {
            "params": {"measurement_and_reporting": {"NAME": "default"}},
            "charges": ["B"],
            "molecules": [],
            "evidence_score_field": None,
        },
        "output": {"charges": []},
    },
]


def adaptor_test():
    for test_dict in TESTS:
        yield adaptor_check, test_dict


def adaptor_check(test_dict):

    result_dict = _parse_evidence_and_format_fixed_labels(test_dict["input"])
    print(result_dict)
    for key in test_dict["output"].keys():
        if key == "fixed_labels":
            assert sorted(test_dict["output"][key].keys()) == sorted(
                result_dict[key].keys()
            )
            for mod_aa in test_dict["output"][key].keys():
                for fl in test_dict["output"][key][mod_aa]:
                    assert fl in result_dict[key][mod_aa]
        elif key in ["charges", "molecules"]:
            assert sorted(test_dict["output"][key]) == sorted(result_dict[key])
        elif key == "params":
            for param_key, param_value in test_dict["output"]["params"].items():
                assert result_dict["params"][param_key] == param_value


if __name__ == "__main__":
    for test_dict in TESTS:
        adaptor_check(test_dict)
