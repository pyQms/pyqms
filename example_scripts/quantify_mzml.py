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
import pyqms
import sys
import pickle
import os

try:
    import pymzml
except:
    print("Please install pymzML via: pip install pymzml")


def main(mzml=None):
    """
    Simple script as template for quantification using pyQms.

    Use e.g. the BSA1.mzML example file. Please download it first using
    'get_example_BSA_file.py'

    Usage:
        ./quantify_mzml.py  mzml_file

    Note:

        The peptides under molecules are BSA peptides.

    """
    molecules = ["HLVDEPQNLIK", "YICDNQDTISSK", "DLGEEHFK"]
    charges = [2, 3, 4, 5]
    metabolic_labels = None
    fixed_labels = None

    lib = pyqms.IsotopologueLibrary(
        molecules=molecules,
        charges=charges,
        metabolic_labels=metabolic_labels,
        fixed_labels=fixed_labels,
        # params           = params,
        verbose=True,
    )
    run = pymzml.run.Reader(mzml)
    mzml_basename = os.path.basename(mzml)
    results = None
    for spectrum in run:
        # print(spectrum.ID)
        scan_time = spectrum.get("MS:1000016")
        if spectrum["ms level"] == 1:
            results = lib.match_all(
                mz_i_list=spectrum.centroidedPeaks,
                file_name=mzml_basename,
                spec_id=spectrum["id"],
                spec_rt=scan_time,
                results=results,
            )
    # pickle.dump(
    #     results,
    #     open(
    #         '{0}_pyQms_results.pkl'.format(mzml_basename),
    #         'wb'
    #     )
    # )
    print(results)
    return


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
    else:
        main(mzml=sys.argv[1])
