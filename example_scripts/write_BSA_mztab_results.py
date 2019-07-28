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

    Warning:

        This example script is specifically for the BSA1.mzML quantification
        results, since file specific meta data is passed. Please use
        'write_mztab_results.py' for a more general script to produce mzTab
        results.

    Note:

        Please note that the ouput in mzTab format is still in beta stage.
        Since pyQms is a raw quantification tool, some meta data has to be
        passed/set manually by the user.

    """
    results_class = pickle.load(open(result_pkl, "rb"))
    # provide meta data as lists of mztab specific formats. Pass directly
    # mztab correct format.

    mztab_meta_info = {
        "protein_search_engine_score": [],
        "psm_search_engine_score": ["[MS,MS:1001475,OMSSA:evalue, ]"],
        "fixed_mod": ["[UNIMOD, UNIMOD:4, Carbamidomethyl, ]"],
        "variable_mod": ["[UNIMOD, UNIMOD:35, Oxidation, ]"],
        "study_variable-description": ["Standard BSA measurement"],
        "ms_run-location": ["BSA1.mzML"],
    }

    results_class.lookup["mztab_meta_info"] = mztab_meta_info

    results_class.write_result_mztab(
        output_file_name="{0}_results.mztab".format(result_pkl)
    )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
    else:
        main(result_pkl=sys.argv[1])
