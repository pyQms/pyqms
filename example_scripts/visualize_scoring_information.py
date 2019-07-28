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
import pprint
from collections import defaultdict as ddict

try:
    import pymzml
except:
    print("Please install pymzML via: pip install pymzml")


def main(mzml=None):
    """
    Example script fort visualizing the m/z and intensity error, which is the
    basis for the scoring of the matches in pyQms.

    Use spectrum 1165 of the BSA1.mzML example file. A subrange of the spectrum
    from m/z 400 to 500 is used.

    Usage:
        ./visualize_scoring_information.py

    Note:
        This example does not require a reader to access MS spectra, since a
        simnple peak list is used.

    """

    peak_list = [
        (404.2492407565097, 2652.905029296875),
        (405.3003310237508, 4831.56103515625),
        (408.8403673369115, 23153.7109375),
        (409.17476109421705, 10182.2822265625),
        (409.5098740355617, 4770.97412109375),
        (411.17196124490727, 3454.364013671875),
        (413.26627826402705, 6861.84912109375),
        (419.3157903165357, 90201.5625),
        (420.2440507067882, 11098.4716796875),
        (420.31917273788645, 22288.9140625),
        (420.73825281590496, 8159.7099609375),
        (421.2406187369968, 3768.656494140625),
        (427.3787652898548, 5680.43212890625),
        (433.3316647490907, 8430.30859375),
        (434.705984428002, 25924.38671875),
        (435.2080179219357, 11041.2060546875),
        (443.6708762397708, 4081.282470703125),
        (443.69049198141124, 5107.13330078125),
        (443.6974813419733, 9135.3125),
        (443.7112735313511, 2517650.0),
        (443.7282222289076, 5571.26025390625),
        (443.7379762316008, 5227.4033203125),
        (444.1998579474954, 3021.341796875),
        (444.21248374593875, 1156173.75),
        (444.71384916266277, 336326.96875),
        (445.21533524843596, 58547.0703125),
        (445.71700965093, 4182.04345703125),
        (446.1200302053469, 93216.3359375),
        (447.09963627699824, 3806.537109375),
        (447.1169242266495, 59846.37109375),
        (447.3464079857604, 13170.9541015625),
        (448.11566395552086, 9294.5107421875),
        (448.3500303628631, 3213.052490234375),
        (452.1123280000919, 5092.0869140625),
        (461.1934526664677, 4022.537353515625),
        (462.1463969367603, 99732.5),
        (463.14561508666384, 24247.015625),
        (464.1433022096936, 20417.041015625),
        (465.1421080732791, 3222.4052734375),
        (470.1669593722212, 8621.81640625),
        (475.23989190282134, 3369.073974609375),
        (493.27465300375036, 2725.885986328125),
        (496.0077303201583, 8604.0830078125),
    ]
    print("{0:-^100}".format("Library generation"))
    lib = pyqms.IsotopologueLibrary(
        molecules=["DDSPDLPK"],
        charges=[2],
        metabolic_labels=None,
        fixed_labels=None,
        verbose=True,
    )
    print("{0:-^100}".format("Library generation"))

    results = lib.match_all(
        mz_i_list=peak_list,
        file_name="BSA_test",
        spec_id=1165,
        spec_rt=29.10,
        results=None,
    )
    for key, i, entry in results.extract_results():
        p = pymzml.plot.Factory()
        label_mz_error = []
        label_i_error = []
        measured_peaks = []
        matched_peaks = []
        peak_info = ddict(list)
        # pprint.pprint(entry.peaks)
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
                mz_error = (measured_mz - calculated_mz) / (measured_mz * 1e-6)
                label_mz_error.append(
                    (calculated_mz, "{0:5.3f} ppm m/z error".format(mz_error))
                )
                scaled_intensity = calculated_intensity * entry.scaling_factor
                rel_i_error = (
                    abs(measured_intensity - scaled_intensity) / scaled_intensity
                )

                peak_info["measured peaks"].append(measured_mz)
                peak_info["theoretical peaks"].append(calculated_mz)
                peak_info["relative intensity"].append(relative_i)
                peak_info["scaled matched peaks"].append(
                    calculated_intensity * entry.scaling_factor
                )
                peak_info["mz error"].append(mz_error)
                peak_info["i error"].append(rel_i_error)

                if rel_i_error > 1:
                    rel_i_error = 1

                label_i_error.append(
                    (calculated_mz, "{0:5.3f} rel. intensity error".format(rel_i_error))
                )

        mz_only = [n[0] for n in measured_peaks]
        mz_range = [min(mz_only) - 1, max(mz_only) + 1]
        peptide = results.lookup["formula to molecule"][key.formula][0]
        p.newPlot(
            header="Formula: {0}; Peptide: {1}; Charge: {2}\n Amount: {3:1.3f}; Score: {4:1.3f}".format(
                key.formula, peptide, key.charge, entry.scaling_factor, entry.score
            ),
            mzRange=mz_range,
        )
        p.add(measured_peaks, color=(0, 0, 0), style="sticks")
        p.add(matched_peaks, color=(0, 200, 0), style="triangles")
        p.add(label_mz_error, color=(255, 0, 0), style="label_x")
        p.add(label_i_error, color=(255, 0, 0), style="label_x")

        plot_name = os.path.join(
            os.pardir,
            "data",
            "Score_visualization_Peptide_{1}_Charge_{2}.xhtml".format(
                key.file_name, peptide, key.charge
            ),
        )
        p.save(filename=plot_name, mzRange=mz_range)
        print("Plotted file {0}".format(plot_name))
        # print(entry)
        print("Match info")
        for key, value_list in sorted(peak_info.items()):
            print(key)
            print("[{0}]".format(",".join([str(n) for n in value_list])))
            print()
    return


if __name__ == "__main__":
    main()
