#!/usr/bin/env python
# encoding: utf-8
"""
Test if the fixed_labeled is counted into the metabolic_labels
given special cases like:

            'fixed_labels' : {'C' : ['C(2) H(3) N(1) O(1)'] },
should give the same results assert
            'fixed_labels'  : {'C' : ['C(2) H(3) 14N(1) O(1)'] },


"""
import pyqms


TESTS = {
    # 2 isotope element (N,nitrogen)
    "14N_is_14N": {
        "set1": {
            "molecules": ["KLEINERTEST"],
            "fixed_labels": {"C": ["C(2) H(3) N(1) O(1)"]},
        },
        "set2": {
            "molecules": ["TESTKLEINER"],
            "fixed_labels": {"C": ["C(2) H(3) 14N(1) O(1)"]},
        },
    },
    "15N_is_15N": {
        "set1": {
            "molecules": ["KLEINERTEST"],
            "fixed_labels": {"C": ["C(2) H(3) N(1) O(1)"]},
            "metabolic_labels": {"15N": [0.994]},
        },
        "set2": {
            "molecules": ["TESTKLEINER"],
            "fixed_labels": {"C": ["C(2) H(3) 15N(1) O(1)"]},
            "metabolic_labels": {"15N": [0.994]},
        },
    },
    "15N_SILAC_is_15N": {
        "set1": {
            "molecules": ["KLEINERTEST"],
            "fixed_labels": {"R": ["C(-6) 13C(6) 15N(-4) 15N(4)"]},
            "metabolic_labels": {"15N": [0.994]},
        },
        "set2": {
            "molecules": ["TESTKLEINER"],
            "fixed_labels": {"R": ["C(-6) 13C(6) N(-4) N(4)"]},
            "metabolic_labels": {"15N": [0.994]},
        },
    },
    ## <><><> ##
}


def create_simple_cam_isotopologue_test():
    for test_id, test_dict in TESTS.items():
        yield create_simple_cam_isotopologue, test_id, test_dict


def create_simple_cam_isotopologue(test_id, test_dict):
    print(
        """
        ol diff
    1868.9315502796
    1869.9285851746
    """
    )
    lib_1 = pyqms.IsotopologueLibrary(charges=[2], verbose=False, **test_dict["set1"])
    lib_2 = pyqms.IsotopologueLibrary(charges=[2], verbose=False, **test_dict["set2"])
    formula_1 = list(lib_1.keys())[0]
    formula_2 = list(lib_2.keys())[0]
    # __oOo__
    for label_percentile in lib_1[formula_1]["env"].keys():
        print(lib_1.lookup["formula to molecule"][formula_1])
        print(lib_2.lookup["formula to molecule"][formula_2])
        for n, mass in enumerate(lib_1[formula_1]["env"][label_percentile]["mass"]):
            print(lib_1[formula_1]["env"][label_percentile]["mass"][n])
            print(lib_2[formula_2]["env"][label_percentile]["mass"][n])
            assert (
                lib_1[formula_1]["env"][label_percentile]["mass"][n]
                - lib_2[formula_2]["env"][label_percentile]["mass"][n]
                < 0.000000001
            )
            break
    # new_distribution = lib._recalc_isotopic_distribution(
    #     test_dict['enriched_element'],
    #     test_dict['target_percentile'],
    #     test_dict['enriched_isotope']
    # )
    # print( new_distribution )
    # print( lib.isotopic_distributions[ test_dict['enriched_element'] ])
    # # percentiles and abundance are the same ...
    # percentiles = set()
    # sum_of_abundances = 0
    # for mass, abundance, isoto_pos in new_distribution:
    #     percentiles.add( uniform( abundance ))
    #     sum_of_abundances += abundance
    # # this uniform function rounds on 2 digits
    # # because _recalc_isotopic_distribution is scaling original
    # # and natural abundances first
    # # errors only occur in lower abundances
    # assert uniform(test_dict['target_percentile']) in percentiles
    # assert 1 - sum_of_abundances <= sys.float_info.epsilon


if __name__ == "__main__":
    for k, v in TESTS.items():
        create_simple_cam_isotopologue(k, v)
