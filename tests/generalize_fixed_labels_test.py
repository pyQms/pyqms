#!/usr/bin/env python3
# encoding: utf-8
'''

Test _extend_kb_with_fixed_labels from core

'''
import pyqms
import sys
import unittest
import copy

TESTS = [  
    # test new generalize fixed labels function which accepts also unimods strings
    {
        'params': {
            'molecules': ['KLEINERTEST'],
            'charges': [2, ],
            'metabolic_labels': {'15N':[0]}
        },
        'fixed_labels_old_style':{
            'R': ['C(-6) 13C(6) N(-4) 15N(4)', ''],
            'K': ['C(-6) 13C(6) N(-2) 15N(2)', '']
        },
        'fixed_labels_new_style':{
            'R': ['Label:13C(6)15N(4)', ''],
            'K': ['Label:13C(6)15N(2)', '']
        },
    },
    # 15N labeling or unlabeled
    {
        'params': {
            'molecules': ['KLEINERCHECK'],
            'charges': [2, ],
            'metabolic_labels': {'15N':[0,0.99]}
        },
        'fixed_labels_old_style':{
            'R': ['C(-6) 13C(6) N(-4) 15N(4)', ''],
            'K': ['C(-6) 13C(6) N(-2) 15N(2)', ''],
            'C': ['O(1) H(3) 14N(1) C(2)'],
        },
        'fixed_labels_new_style':{
            'R': ['Label:13C(6)15N(4)', ''],
            'K': ['Label:13C(6)15N(2)', ''],
            'C': ['Carbamidomethyl']
        },
    },
    {
        'params': {
            'molecules': ['KLEINERCHECK'],
            'charges': [2, ],
        },
        'fixed_labels_old_style':{
            'R': ['C(-6) 13C(6) N(-4) 15N(4)', ''],
            'K': ['C(-6) 13C(6) N(-2) 15N(2)', ''],
            'C': ['O(1) H(3) N(1) C(2)'],
        },
        'fixed_labels_new_style':{
            'R': ['Label:13C(6)15N(4)', ''],
            'K': ['Label:13C(6)15N(2)', ''],
            'C': ['Carbamidomethyl']
        },
    },
    # 13C labeling or unlabeled
    {
        'params': {
            'molecules': ['KLEINERCHECK'],
            'charges': [2, ],
            'metabolic_labels': {'13C':[0,0.99]}
        },
        'fixed_labels_old_style':{
            'R': ['C(-6) 13C(6) N(-4) 15N(4)', ''],
            'K': ['C(-6) 13C(6) N(-2) 15N(2)', ''],
            'C': ['O(1) H(3) N(1) 12C(2)'],
        },
        'fixed_labels_new_style':{
            'R': ['Label:13C(6)15N(4)', ''],
            'K': ['Label:13C(6)15N(2)', ''],
            'C': ['Carbamidomethyl']
        },
    },
    {
        'params': {
            'molecules': ['KLEINERCHECK'],
            'charges': [2, ],
            'metabolic_labels': {'13C':[0]}
        },
        'fixed_labels_old_style':{
            'R': ['C(-6) 13C(6) N(-4) 15N(4)', ''],
            'K': ['C(-6) 13C(6) N(-2) 15N(2)', ''],
            'C': ['O(1) H(3) N(1) C(2)'],
        },
        'fixed_labels_new_style':{
            'R': ['Label:13C(6)15N(4)', ''],
            'K': ['Label:13C(6)15N(2)', ''],
            'C': ['Carbamidomethyl']
        },
    },
]


def extend_kb_with_fixed_labels_test():
    for test_id, test_dict in enumerate(TESTS):
        yield _extend_kb_with_fixed_labels, test_id, test_dict


def _extend_kb_with_fixed_labels( test_id, test_dict ):

    params_1= copy.deepcopy(test_dict['params'])
    params_1['fixed_labels'] = test_dict['fixed_labels_old_style']
    lib_1 = pyqms.IsotopologueLibrary(
        **params_1
    )

    params_2= copy.deepcopy(test_dict['params'])
    params_2['fixed_labels'] = test_dict['fixed_labels_new_style']

    lib_2 = pyqms.IsotopologueLibrary(
        **params_2
    )
    assert sorted(lib_1.keys()) == sorted(lib_2.keys())


if __name__ == '__main__':
    pass
