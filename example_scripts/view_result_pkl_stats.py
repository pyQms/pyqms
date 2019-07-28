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


def main(result_pkl=None):
    """

    usage:
        ./view_result_pkl_stats.py <Path2ResultPkl>

    This script will show the stats of a result pkl file. Can be used to query
    the number of quantified formulas and charge states etc.

    """
    results_class = pickle.load(open(result_pkl, "rb"))
    print("Result pkl file holds the following information:")
    print()
    for key, value in results_class.index.items():
        print("Number of {0: <20}: {1}".format(key, len(value)))
        print("\tExample values (up to 5): {0}".format(list(value)[:5]))
        print()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
    else:
        main(result_pkl=sys.argv[1])
