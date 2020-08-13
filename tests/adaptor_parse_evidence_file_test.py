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

from pyqms.adaptors import parse_evidence as parse_evidence
import os
import pyqms

tmp_cc_factory = pyqms.chemical_composition.ChemicalComposition()

tmp_cc_factory.add_chemical_formula("C2O1H3N1")
TESTS = [
    {
        # evidence file parse
        "input": {
            "fixed_labels": None,
            "molecules": None,
            "evidence_score_field": None,
            "evidence_files": [os.path.join("tests", "data", "test_BSA_evidence.csv")],
        },
        "output": {
            "molecules": [
                "CCTESLVNR#Carbamidomethyl:1;Carbamidomethyl:2",
                "CCTESLVNR",
                "DDSPDLPK",
                "+C37H59N9O16",
            ]
        },
    },
    {
        # evidence file parse
        "input": {
            "fixed_labels": {
                "C": [
                    {
                        "element_composition": {"O": 1, "H": 3, "14N": 1, "C": 2},
                        "evidence_mod_name": "Carbamidomethyl",
                    }
                ]
            },
            "molecules": None,
            "evidence_score_field": None,
            "evidence_files": [os.path.join("tests", "data", "test_BSA_evidence.csv")],
        },
        "output": {"molecules": ["CCTESLVNR", "DDSPDLPK", "+C37H59N9O16"]},
    },
    {
        # evidence file parse
        "input": {
            "fixed_labels": {
                "C": [
                    {
                        "element_composition": tmp_cc_factory,
                        "evidence_mod_name": "Carbamidomethyl",
                    }
                ]
            },
            "molecules": None,
            "evidence_score_field": None,
            "evidence_files": [os.path.join("tests", "data", "test_BSA_evidence.csv")],
        },
        "output": {"molecules": ["CCTESLVNR", "DDSPDLPK", "+C37H59N9O16"]},
    },
    {
        # full evidence parse
        "input": {
            "fixed_labels": {
                "C": [
                    {
                        "element_composition": tmp_cc_factory,
                        "evidence_mod_name": "Carbamidomethyl",
                    }
                ]
            },
            "molecules": None,
            "evidence_score_field": None,
            "evidence_files": [os.path.join("data", "BSA1_omssa_2_1_9_unified.csv")],
        },
        "output": {
            "molecules": [
                "AEFVEVTK",
                "CCTESLVNR",
                "DDSPDLPK",
                "DLGEEHFK",
                "EACFAVEGPK",
                "ECCDKPLLEK",
                "ETYGDMADCCEK",
                "EYEATLEECCAK",
                "GACLLPK",
                "HLVDEPQNLIK",
                "LCVLHEK",
                "LKPDPNTLCDEFK",
                "LKPDPNTLCDEFK#Oxidation:1",
                "LKPDPNTLCDEFK#Oxidation:1;Oxidation:2",
                "LVTDLTK",
                "LVVSTQTALA",
                "SHCIAEVEK",
                "YICDNQDTISSK",
                "YLYEIAR",
            ]
        },
    },
    {
        # full evidence parse
        "input": {
            "fixed_labels": {
                "C": [
                    {
                        "element_composition": tmp_cc_factory,
                        "evidence_mod_name": "Carbamidomethyl",
                    }
                ]
            },
            "molecules": None,
            "evidence_score_field": None,
            "evidence_files": [os.path.join("data", "BSA1_omssa_2_1_9.mztab")],
        },
        "output": {
            "molecules": [
                "AEFVEVTK",
                "CCTESLVNR",
                "DDSPDLPK",
                "DLGEEHFK",
                "EACFAVEGPK",
                "ECCDKPLLEK",
                "ETYGDMADCCEK",
                "EYEATLEECCAK",
                "GACLLPK",
                "HLVDEPQNLIK",
                "LCVLHEK",
                "LKPDPNTLCDEFK",
                "LKPDPNTLCDEFK#Oxidation:1",
                "LKPDPNTLCDEFK#Oxidation:1;Oxidation:2",
                "LVTDLTK",
                "LVVSTQTALA",
                "SHCIAEVEK",
                "YICDNQDTISSK",
                "YLYEIAR",
            ]
        },
    },
]


def adaptor_test():
    for test_dict in TESTS:
        yield adaptor_check, test_dict


def adaptor_check(test_dict):

    formatted_fixed_labels, evidence_lookup, molecule_list = parse_evidence(
        **test_dict["input"]
    )
    assert len(molecule_list) == len(test_dict["output"]["molecules"])
    assert sorted(molecule_list) == sorted(test_dict["output"]["molecules"])
    for formula in evidence_lookup.keys():
        for molecule in evidence_lookup[formula].keys():
            if len(evidence_lookup[formula][molecule]['trivial_names']) !=0:
                for ev_dict in evidence_lookup[formula][molecule]['evidences']:
                    assert 'trivial_name' in ev_dict.keys()
                    trivial_name = ev_dict['trivial_name']
                    if ';' in trivial_name:
                        trivial_name_list = trivial_name.split(';')
                    else:
                        trivial_name_list = [ trivial_name ]
                    for tmp_trivial_name in trivial_name_list:
                        assert tmp_trivial_name in evidence_lookup[formula][molecule]['trivial_names']
                    assert sorted(evidence_lookup[formula][molecule]['trivial_names']) == sorted(trivial_name_list)

if __name__ == "__main__":
    for test_dict in TESTS:
        adaptor_check(test_dict)
