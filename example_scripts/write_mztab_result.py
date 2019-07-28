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
        ./write_mztab_results.py <Path2ResultPkl>

    Will write all results of a result pkl into a .mztab file. Please refer to
    Documentation of :doc:`results` for further information.

    Note:

        Please note that the ouput in mzTab format is still in beta stage.
        Since pyQms is a raw quantification tool, some meta data has to be
        passed/set manually by the user.

    """
    results_class = pickle.load(open(result_pkl, "rb"))

    results_class.write_result_mztab(
        output_file_name="{0}_results.mztab".format(result_pkl)
    )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
    else:
        main(result_pkl=sys.argv[1])
