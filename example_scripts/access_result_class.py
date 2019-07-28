#!/usr/bin/env python3
# encoding: utf-8
"""
    pyQms
    -----

    Python module for fast and accurate mass spectrometry data quantification

    :license: MIT, see LICENSE.txt for more details

    Authors:

        * Leufken, J.
        * Niehues, A.
        * Sarin, L.P.
        * Hippler, M.
        * Leidel, S.A.
        * Fufezan, C.

"""
import pickle
import sys
import pprint


def main(result_pkl=None):
    """

    usage:
        ./access_result_class.py <Path2ResultPkl>

    This script will produce a dictionary with all summed up peptide amounts.
    The main idea is to show how to access the result pkl and loop over the
    data structure.

    Note:

        Since no filters (score, RT windows, etc.) are applied, this script
        should not be used to estimate the actual amount of the quantified
        molecules in the results pkl.


    """
    results_class = pickle.load(open(result_pkl, "rb"))

    amount_collector = {}

    for key, value in results_class.items():
        peptide = results_class.lookup["formula to molecule"][key.formula][0]
        if peptide not in amount_collector.keys():
            amount_collector[peptide] = {"amount": 0}
        for matched_spectrum in value["data"]:
            amount_collector[peptide]["amount"] += matched_spectrum.scaling_factor

    pprint.pprint(amount_collector)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
    else:
        main(result_pkl=sys.argv[1])
