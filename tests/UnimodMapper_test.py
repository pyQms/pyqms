#!/usr/bin/env python
# encoding: utf-8

import pyqms
import unittest

M = pyqms.UnimodMapper()

UNIMODMAPPER_FUNCTIONS = [
    M.name2mass,
    M.name2composition,
    M.name2id,
    M.id2mass,
    M.id2composition,
    M.id2name,
    M.mass2name_list,
    M.mass2id_list,
    M.mass2composition_list,
    M.appMass2id_list,
    M.appMass2element_list,
    M.appMass2name_list,
    M.composition2name_list,
    M.name2specificity_site_list,
    M.composition2id_list,
    M.composition2mass,
    M._map_key_2_index_2_value,
]

TESTS = [
    [
        {
            'in' : {
                'args': ['ICAT-G:2H(8)']
            },
            'out': 494.30142
        },  # First test : pyqms.UnimodMapper.name2mass,
        {
            'in' : {
                'args': ['ICAT-G:2H(8)']
            },
            'out': {'H': 30, 'C': 22, 'O': 6, 'S': 1, '2H': 8, 'N': 4}
        },  # pyqms.UnimodMapper.name2composition,
        {
            'in' : {
                'args': ['ICAT-G:2H(8)']
            },
            'out' : '9'
        },  # pyqms.UnimodMapper.name2id,
        {
            'in' : {
                'args': ['9']
            },
            'out': 494.30142
        },  # pyqms.UnimodMapper.id2mass,
        {
            'in' : {
                'args': ['9']
            },
            'out' : {'N': 4, 'S': 1, '2H': 8, 'O': 6, 'C': 22, 'H': 30}
        },  # pyqms.UnimodMapper.id2composition,
        {
            'in' : {
                'args': '9'
            },
            'out' : 'ICAT-G:2H(8)'
        },  # pyqms.UnimodMapper.id2name,
        {
            'in' : {
                'args': [494.30142]
            },
            'out': ['ICAT-G:2H(8)']
        },  # pyqms.UnimodMapper.mass2name_list,
        {
            'in' : {
                'args': [494.30142]
            },
            'out' : ['9']
        },  # pyqms.UnimodMapper.mass2id_list,
        {
            'in' : {
                'args': [494.30142]
            },
            'out' : [{'N': 4, 'S': 1, '2H': 8, 'O': 6, 'C': 22, 'H': 30}]
        },  # pyqms.UnimodMapper.mass2composition_list,
        {
            'in' : {
                'args': [18],
                'kwargs': {'decimal_places': 0 }
            },
            'out': ['127', '329', '608', '1079', '1167']
        },  # pyqms.UnimodMapper.appMass2id_list
        {
            'in' : {
                'args': [18],
                'kwargs': {'decimal_places': 0 }
            },
            'out': [
                {'F': 1, 'H': -1}, {'13C': 1, 'H': -1, '2H': 3},
                {'H': -2, 'C': -1, 'S': 1},
                {'H': 2, 'C': 4, 'O': -2},
                {'H': -2, 'C': -1, 'O': 2}
            ]
        },  # pyqms.UnimodMapper.appMass2element_list
        {
            'in' : {
                'args': [18],
                'kwargs': {'decimal_places': 0 }
            },
            'out': [
                'Fluoro',
                'Methyl:2H(3)13C(1)',
                'Xle->Met',
                'Glu->Phe',
                'Pro->Asp'
            ]
        },  # pyqms.UnimodMapper.appMass2name_list
        {
            'in' : {
                'args': ['C(2)H(3)N(1)O(1)'],
            },
            'out': ['Carbamidomethyl', 'Ala->Gln', 'Gly->Asn', 'Gly']
        },
        {
            'in' : {
                'args': ['Ala->Gln'],
            },
            'out': ['A']
        },
        {
            'in' : {
                'args': [ 'C(22)H(30)2H(8)N(4)O(6)S(1)'],
            },
            'out': ['9']
        },
        {
            'in' : {
                'args': [ 'C(22)H(30)2H(8)N(4)O(6)S(1)'],
            },
            'out': 494.30142
        },
        {
            'in' : {
                'args': [ 'ThisKeyIsNotPresent', 'mass'],
            },
            'out': None
        },

        #_map_key_2_index_2_value


        # end of first data set ...
    ]
]


def test_set_integirty_test():
    for test_id, list_of_test in enumerate(TESTS):
        yield input_list_check, list_of_test


def input_list_check( list_of_test ):
    '''
    Checks that the number of tests can actually be zipped
    with the unimodMapper functions
    '''
    assert len(list_of_test) <= len(UNIMODMAPPER_FUNCTIONS)


def unimodMapper_conversion_test():
    for test_id, list_of_test in enumerate(TESTS):
        test_function_association = zip( list_of_test, UNIMODMAPPER_FUNCTIONS )
        for test_dict, mapper_function in test_function_association:
            if 'kwargs' not in test_dict['in'].keys():
                test_dict['in']['kwargs'] = {}
            yield mapper_function_check, test_dict, mapper_function


def mapper_function_check( test_dict, mapper_function):
        mapper_output = mapper_function(
            *test_dict['in']['args'],
            **test_dict['in']['kwargs']
        )
        print(mapper_output, test_dict['out'])
        assert mapper_output == test_dict['out']



class TestResults( unittest.TestCase ):
    def setUp( self ):
        self.alt_mapper = pyqms.UnimodMapper()
        self.alt_mapper.unimod_xml_name = 'wrong_unimod_xml_name.xml'

    def crash_test(self):
        with self.assertRaises(SystemExit) as system_exit_check:
            self.alt_mapper._parseXML()
        self.assertEqual(
            system_exit_check.exception.code,
            1
        )
        return


if __name__ == '__main__':
    print('Yes!')
