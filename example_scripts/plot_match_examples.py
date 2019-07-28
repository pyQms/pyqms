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
    import pymzml
    import pymzml.plot
except:
    print("Please install pymzML via: pip install pymzml")


def main(result_pkl=None):
    """

    usage:
        ./plot_match_examples.py <Path2ResultPkl>

    Extracts the match information and plots one example isotopologue match into
    the 'data' folder. Uses the plot function of pymzML (`pymzML.plot`_). Use
    this script as template for annotating spectra with match information.

    Note:

        Plots only one high scored formula (mScore >0.95) from the result pkl.
        Use e.g. with the 'BSA1.mzML_pyQms_results.pkl' obtained from e.g.
        example script parse_ident_file_and_quantify_with_carbamidomethylation.py
        to get example plotting data.

    .. _pymzML.plot:
        https://pymzml.github.io/plot.html

    """
    results_class = pickle.load(open(result_pkl, "rb"))

    for key, i, entry in results_class.extract_results():
        if entry.score > 0.95:
            p = pymzml.plot.Factory()
            label_x = []
            measured_peaks = []
            matched_peaks = []
            for (
                measured_mz,
                measured_intensity,
                relative_i,
                calculated_mz,
                calculated_intensity,
            ) in entry.peaks:
                if measured_mz is not None:
                    measured_peaks.append((measured_mz, measured_intensity))
                    matched_peaks.append(
                        (calculated_mz, calculated_intensity * entry.scaling_factor)
                    )
                    label_x.append(
                        (
                            calculated_mz,
                            "{0:5.3f} ppm".format(
                                (measured_mz - calculated_mz) / (measured_mz * 1e-6)
                            ),
                        )
                    )

            mz_only = [n[0] for n in measured_peaks]
            mz_range = [min(mz_only) - 1, max(mz_only) + 1]
            peptides = results_class.lookup["formula to molecule"][key.formula]
            if len(peptides) > 1:
                continue
            p.newPlot(
                header="Formula: {0}; Peptide: {1}; Charge: {2}\n File: {3}; Scan: {4}; RT: {5:1.3f}\n Amount: {6:1.3f}; Score: {7:1.3f}".format(
                    key.formula,
                    peptides[0],
                    key.charge,
                    key.file_name,
                    entry.spec_id,
                    entry.rt,
                    entry.scaling_factor,
                    entry.score,
                ),
                mzRange=mz_range,
            )
            p.add(measured_peaks, color=(0, 0, 0), style="sticks")
            p.add(matched_peaks, color=(0, 200, 0), style="triangles")
            p.add(label_x, color=(0, 0, 255), style="label_x")

            plot_name = os.path.join(
                os.pardir,
                "data",
                "{0}_Peptide_{1}_Charge_{2}.xhtml".format(
                    key.file_name, peptides[0], key.charge
                ),
            )
            p.save(filename=plot_name, mzRange=mz_range)
            print("Plotted file {0}".format(plot_name))
            break


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(main.__doc__)
    else:
        main(result_pkl=sys.argv[1])
