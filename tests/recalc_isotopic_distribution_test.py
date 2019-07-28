#!/usr/bin/env python
# encoding: utf-8
"""

Test _recalc_isotopic_distribution from core



"""
import pyqms
import sys


TESTS = {
    # 2 isotope element (N,nitrogen)
    "mo_test1": {
        "enriched_element": "N",
        "enriched_isotope": 15,
        "target_percentile": 0.994,
    },
    # 3 isotope element (oxygen)
    "mo_test2": {
        "enriched_element": "O",
        "enriched_isotope": 18,
        "target_percentile": 0.970,
    },
    "mo_test3": {
        "enriched_element": "N",
        "enriched_isotope": 15,
        "target_percentile": 0.32,
    },
    # 'mo_test4' : {
    #     'enriched_element'  : 'N',
    #     'enriched_isotope'  : 14,
    #     'target_percentile' : 0.994,
    # }
}


def uniform(abundance):
    return str(round(abundance, 3))


def recalc_isotopic_distribution_test():
    for test_id, test_dict in TESTS.items():
        yield _recalc_isotopic_distribution, test_id, test_dict


def _recalc_isotopic_distribution(test_id, test_dict):
    lib = pyqms.IsotopologueLibrary(
        molecules=["KLEINERTEST"], charges=[2], verbose=False
    )
    new_distribution = lib._recalc_isotopic_distribution(
        element=test_dict["enriched_element"],
        target_percentile=test_dict["target_percentile"],
        enriched_isotope=test_dict["enriched_isotope"],
    )
    print(new_distribution)
    print(lib.isotopic_distributions[test_dict["enriched_element"]])
    # percentiles and abundance are the same ...
    percentiles = set()
    sum_of_abundances = 0
    for mass, abundance, isoto_pos in new_distribution:
        percentiles.add(uniform(abundance))
        sum_of_abundances += abundance
    # this uniform function rounds on 2 digits
    # because _recalc_isotopic_distribution is scaling original
    # and natural abundances first
    # errors only occur in lower abundances
    assert uniform(test_dict["target_percentile"]) in percentiles
    assert 1 - sum_of_abundances <= sys.float_info.epsilon


if __name__ == "__main__":
    for k, v in TESTS.items():
        _recalc_isotopic_distribution(k, v)
