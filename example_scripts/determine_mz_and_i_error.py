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
import os


def main(result_pkl=None):
    """

    usage:
        ./determine_mz_and_i_error.py <Path2ResultPkl>

    This script will determine the apparant mz and intensity error present
    in the quantifications for the given result pkl.


    """
    results_class = pickle.load(open(result_pkl, "rb"))

    plot_name = os.path.join(
        os.path.dirname(result_pkl),
        "mz_and_intensity_error_{0}.pdf".format(os.path.basename(result_pkl)),
    )

    results_class._determine_measured_error(
        score_threshold=None, topX=3, filename=plot_name, plot=True
    )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
    else:
        main(result_pkl=sys.argv[1])
