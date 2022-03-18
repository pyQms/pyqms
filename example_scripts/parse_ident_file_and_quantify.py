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
import pyqms.adaptors

try:
    import pymzml
except:
    print("Please install pymzML via: pip install pymzml")


def main(ident_file=None, mzml_file=None):
    """
    Script to automatically parse `Ursgal`_ result files and quantify it via
    pyQms. Please refer to Documenation of :doc:`adaptors` for further
    information.

    `Ursgal`_ result files or files in `mzTab` format are read in and used for
    quantification of the BSA example file.

    Note:

        Use e.g. the BSA1.mzML example file. Please download it first using
        'get_example_BSA_file.py'. Evidence files can also be found in the
        data folder 'BSA1_omssa_2_1_9_unified.csv' or 'BSA1_omssa_2_1_9.mztab'

    Usage:

        ./parse_ident_file_and_quantify.py <ident_file> <mzml_file>

    .. _Ursgal:
        https://github.com/ursgal/ursgal

    .. _mzTab:
        http://www.psidev.info/mztab

    """

    if ident_file.upper().endswith("MZTAB"):
        evidence_score_field = "search_engine_score[1]"
    else:
        # this is the default value in the adaptor
        evidence_score_field = "PEP"

    print('Evidence score field "{0}" will be used.'.format(evidence_score_field))

    fixed_labels, evidences, molecules = pyqms.adaptors.parse_evidence(
        fixed_labels=None,
        evidence_files=[ident_file],
        evidence_score_field=evidence_score_field,
    )

    params = {
        "molecules": molecules,
        "charges": [1, 2, 3, 4, 5],
        "metabolic_labels": {"15N": [0]},
        "fixed_labels": fixed_labels,
        "verbose": True,
        "evidences": evidences,
    }

    lib = pyqms.IsotopologueLibrary(**params)

    run = pymzml.run.Reader(mzml_file)
    out_folder = os.path.dirname(mzml_file)
    mzml_file_basename = os.path.basename(mzml_file)
    results = None
    for spectrum in run:
        try:
            # pymzML 2.0.0 style
            scan_time = spectrum.scan_time
        except:
            # scan time will be in seconds
            scan_time = spectrum.get("MS:1000016")
        if spectrum["ms level"] == 1:
            results = lib.match_all(
                mz_i_list=spectrum.centroidedPeaks,
                file_name=mzml_file_basename,
                spec_id=spectrum["id"],
                spec_rt=scan_time,
                results=results,
            )
    pickle.dump(
        results,
        open(
            os.path.join(
                out_folder, "{0}_pyQms_results.pkl".format(mzml_file_basename)
            ),
            "wb",
        ),
    )
    return


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
    else:
        main(ident_file=sys.argv[1], mzml_file=sys.argv[2])
