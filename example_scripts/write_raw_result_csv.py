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
        ./write_raw_result_csv.py <Path2ResultPkl>

    Will write all results of a result pkl into a .csv file. Please refer to
    Documentation of :doc:`results` for further information.

    Warning:

        The resulting .csv files can become very large depending on the provided
        pkl file!

    Keys in csv:

        * Formula           : molecular formula of the molecule (str)
        * Molecule          : molecule or trivial name (str)
        * Charge            : charge of the molecule (int)
        * ScanID            : ScanID of the quantified spectrum (int)
        * Label Percentiles : Labeling percentile ( (element, enrichment in %), )
        * Amount            : the determined amount of the molecule
        * Retention Time    : retetention time of the ScanID
        * mScore            : score of the isotopologue match
        * Filename          : filename of spectrum input files

    """
    results_class = pickle.load(open(result_pkl, "rb"))

    results_class.write_result_csv(
        output_file_name="{0}_raw_results.csv".format(result_pkl)
    )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
    else:
        main(result_pkl=sys.argv[1])
