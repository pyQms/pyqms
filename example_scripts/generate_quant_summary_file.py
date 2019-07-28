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
        ./generate_quant_summary_file.py <Path2ResultPkl>

    This script will produce quant summary file with all according evidence
    information which are stored in the result pkl file. Amounts (maxI) will be
    calculated if possible.

    Note:

        Make sure, an evidence lookup is provided in the results class, so that
        retention time windows can be defined. Otherwise no meaningful amounts
        can be calculated.


    Warning:

        Can take very long depending on pkl size!

    """
    results_class = pickle.load(open(result_pkl, "rb"))
    rt_border_tolerance = 1
    # quant_summary_file  = '{0}_quant_summary.csv'.format(result_pkl)
    quant_summary_file = "{0}_quant_summary.xlsx".format(result_pkl)
    results_class.write_rt_info_file(
        output_file=quant_summary_file,
        list_of_csvdicts=None,
        trivial_name_lookup=None,
        rt_border_tolerance=rt_border_tolerance,
        update=True,
    )
    results_class.calc_amounts_from_rt_info_file(
        rt_info_file=quant_summary_file,
        rt_border_tolerance=rt_border_tolerance,
        calc_amount_function=None,
    )
    return


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
    else:
        main(result_pkl=sys.argv[1])
