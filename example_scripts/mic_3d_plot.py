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

try:
    import rpy2
except:
    print(
        "rpy2 is not installed but required for plotting, please install it and try again"
    )
    print("pip3.4 install rpy2")


def main(pickle_file):
    """
    usage:
        ./mic_3d_plot.py <path_to_pickled_result_class>

    Plots 3-dimensional matched isotope chromatograms (MICs) of pyQms
    quantification results.

    Pickled result class can contain thousands of molecules therefore this
    example script stops plotting after 10 plotted MICs. Otherwise all
    quantified formula-charge-filename combinations will be plotted!

    Use e.g. the BSA data example. Download via 'get_example_BSA_file.py' and
    quantify using 'parse_ident_file_and_quantify_with_carbmidomethylation.py'.

    Note:

        Installation of R and rpy2 is required.

    """
    results = pickle.load(open(pickle_file, "rb"))
    out_folder = os.path.join(os.path.dirname(pickle_file), "plots")
    if os.path.exists(out_folder) is False:
        os.mkdir(out_folder)
    print("Plotting into folder: {0}".format(out_folder))
    if len(results.keys()) > 10:
        print(
            """
Result class should not hold more then 10 keys, to prevent plot overflow!
Will stop after 10 plots!
"""
        )
        # sys.exit()
    for n, key in enumerate(results.keys()):
        if n > 10:
            print("Stopping after 10 plots!")
            exit()
        if len(results[key]["data"]) <= 15:
            continue
        mzml_filename = key.file_name
        if os.sep in mzml_filename:
            mzml_filename = os.path.basename(mzml_filename)

        file_name = os.path.join(
            out_folder,
            "MIC_3D_{0}_{1}_{2}_{3}".format(
                "_".join(results.lookup["formula to molecule"][key.formula]),
                key.charge,
                key.label_percentiles,
                mzml_filename,
            ),
        )
        results.plot_MIC_3D(key, file_name=file_name)

    return


if __name__ == "__main__":
    if len(sys.argv) <= 1:
        sys.exit(main.__doc__)
    main(sys.argv[1])
