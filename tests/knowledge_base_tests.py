#!/usr/bin/env python
# encoding: utf-8
"""

Testfunction to test the integrity of the knowledge_base.py

Implements some basic tests to check the data of knowledge_base.py

"""
import sys
import pyqms


def element_isotopic_abundance_test():
    for (
        element,
        isotope_atomic_mass_abundance_list,
    ) in pyqms.knowledge_base.isotopic_distributions.items():
        yield check_abundance_sum, isotope_atomic_mass_abundance_list


def check_abundance_sum(isotopic_distribution_list):
    sum_of_abundances = sum([abun for mass, abun in isotopic_distribution_list])
    assert 1 - sum_of_abundances <= sys.float_info.epsilon


if __name__ == "__main__":
    element_isotopic_abundance_test()
