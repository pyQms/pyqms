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
    pyQms.

    For evidence files with molecules with Caramidomethylation as fixed
    modification. These mode will be stripped from the molecules. This is
    important if an metabolic label (like 15N) is applied. This ensures that the
    nitrogens pools of the peptides (which are 15N labeled) do not mix up with
    the nitrogen pool of the Carbamidomethylation (14N since intriduced during
    sample preparation). Please refer to Documenation of :doc:`adaptors` for
    further information.

    `Ursgal`_ result files or files in `mzTab` format are read in and used for
    quantification of the BSA example file.

    Note:

        Use e.g. the BSA1.mzML example file. Please download it first using
        'get_example_BSA_file.py'. Evidence files can also be found in the
        data folder 'BSA1_omssa_2_1_9_unified.csv' or 'BSA1_omssa_2_1_9.mztab'

    Usage:

        ./parse_ident_file_and_quantify_with_carbamidomethylation.py <ident_file> <mzml_file>

    .. _Ursgal:
        https://github.com/ursgal/ursgal

    .. _mzTab:
        http://www.psidev.info/mztab

    """

    # define the fixed label for Caramidomethyl
    tmp_fixed_labels = {
        "C": [
            {
                "element_composition": {"O": 1, "H": 3, "14N": 1, "C": 2},
                "evidence_mod_name": "Carbamidomethyl",
            }
        ]
    }

    formatted_fixed_labels, evidence_lookup, molecule_list = pyqms.adaptors.parse_evidence(
        fixed_labels=tmp_fixed_labels, evidence_files=[ident_file]
    )

    params = {
        "molecules": molecule_list,
        "charges": [1, 2, 3, 4, 5],
        "metabolic_labels": {"15N": [0]},
        "fixed_labels": formatted_fixed_labels,
        "verbose": True,
        "evidences": evidence_lookup,
    }

    lib = pyqms.IsotopologueLibrary(**params)

    run = pymzml.run.Reader(mzml_file)
    out_folder = os.path.dirname(mzml_file)
    mzml_file_basename = os.path.basename(mzml_file)
    results = None
    for spectrum in run:
        spec_id = spectrum["id"]
        try:
            # pymzML 2.0.0 style
            scan_time, unit = spectrum.scan_time
            if "unit" == "minute":
                scan_time /= 60.0
        except:
            # scan time will be in seconds
            scan_time = spectrum.get("MS:1000016") / 60.0
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
