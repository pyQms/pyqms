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
    Examples script to demonstrate a (example) workflow from mzML files to
    peptide abundances. Will plot for every quantified peptide a matched
    isotopologue chromatogram (MIC). The plots include RT windows, maximum
    amount in RT window and identification RT(s).

    `Ursgal`_ result files or files in `mzTab` format are read in and used for
    quantification of the BSA example file.

    Note:

        Use e.g. the BSA1.mzML example file. Please download it first using
        'get_example_BSA_file.py'. Evidence files can also be found in the
        data folder 'BSA1_omssa_2_1_9_unified.csv' or 'BSA1_omssa_2_1_9.mztab'

    Usage:

        ./complete_BSA_quantification.py <ident_file> <mzml_file>

    .. _Ursgal:
        https://github.com/ursgal/ursgal

    .. _mzTab:
        http://www.psidev.info/mztab

    Note:
        rpy2 is required for all plotting

    """

    # define the fixed label for Carbamidomethyl
    tmp_fixed_labels = {
        "C": [
            {
                "element_composition": {"O": 1, "H": 3, "14N": 1, "C": 2},
                "evidence_mod_name": "Carbamidomethyl",
            }
        ]
    }
    if ident_file.upper().endswith("MZTAB"):
        evidence_score_field = "search_engine_score[1]"
    else:
        # this is the default value in the adaptor
        evidence_score_field = "PEP"

    print('Evidence score field "{0}" will be used.'.format(evidence_score_field))
    formatted_fixed_labels, evidence_lookup, molecule_list = pyqms.adaptors.parse_evidence(
        fixed_labels=tmp_fixed_labels,
        evidence_files=[ident_file],
        evidence_score_field=evidence_score_field,
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
    # print(results)
    out_folder = os.path.join(
        os.path.dirname(ident_file), "complete_BSA_quantification"
    )
    if os.path.exists(out_folder) is False:
        os.mkdir(out_folder)
    print()
    print("All results go into folder: {0}".format(out_folder))
    rt_border_tolerance = 1
    quant_summary_file = os.path.join(
        out_folder, "complete_BSA_quantification_summary.xlsx"
    )
    results.write_rt_info_file(
        output_file=quant_summary_file,
        list_of_csvdicts=None,
        trivial_name_lookup=None,
        rt_border_tolerance=rt_border_tolerance,
        update=True,
    )
    calculated_amounts = results.calc_amounts_from_rt_info_file(
        rt_info_file=quant_summary_file,
        rt_border_tolerance=rt_border_tolerance,
        calc_amount_function=None,  # calc_amount_function
    )
    # print(calculated_amounts)
    formula_charge_to_quant_info = {}
    for line_dict in calculated_amounts:
        formula_charge_to_quant_info[
            (line_dict["formula"], int(line_dict["charge"]))
        ] = {
            "rt": line_dict["max I in window (rt)"],
            "amount": line_dict["max I in window"],
            "rt start": line_dict["start (min)"],
            "rt stop": line_dict["stop (min)"],
            "evidence_rts": [],
        }
        if (
            len(
                formula_charge_to_quant_info[
                    (line_dict["formula"], int(line_dict["charge"]))
                ]["evidence_rts"]
            )
            == 0
        ):
            for ev_string in line_dict["evidences (min)"].split(";"):
                formula_charge_to_quant_info[
                    (line_dict["formula"], int(line_dict["charge"]))
                ]["evidence_rts"].append(round(float(ev_string.split("@")[1]), 2))
    import_ok = False
    try:
        import rpy2

        import_ok = True
    except:
        pass
    if import_ok:
        print(
            "Plotting results plot including RT windows, abundances and identifications"
        )
        for key in results.keys():
            short_key = (key.formula, key.charge)

            match_list = results[key]["data"]
            if len(match_list) < 15:
                continue
            file_name = os.path.join(
                out_folder,
                "MIC_2D_{0}_{1}.pdf".format(
                    "_".join(results.lookup["formula to molecule"][key.formula]),
                    key.charge,
                ),
            )
            graphics, grdevices = results.init_r_plot(file_name)

            ablines = {
                key: [
                    {"v": formula_charge_to_quant_info[short_key]["rt"], "lty": 2},
                    {
                        "v": formula_charge_to_quant_info[short_key]["rt start"],
                        "lty": 2,
                        "col": "blue",
                    },
                    {
                        "v": formula_charge_to_quant_info[short_key]["rt stop"],
                        "lty": 2,
                        "col": "blue",
                    },
                ]
            }
            # print(formula_charge_to_quant_info[short_key])
            additional_legends = {
                key: [
                    {
                        "x": formula_charge_to_quant_info[short_key]["rt"],
                        "y": formula_charge_to_quant_info[short_key]["amount"],
                        "text": "max intensity: {0:1.3e}".format(
                            formula_charge_to_quant_info[short_key]["amount"]
                        ),
                        "pos": 3,  # above
                    },
                    {
                        "x": formula_charge_to_quant_info[short_key]["rt start"],
                        "y": formula_charge_to_quant_info[short_key]["amount"] / 2,
                        "text": "RT Window start",
                        "pos": 4,  # right
                        "col": "blue",
                    },
                    {
                        "x": formula_charge_to_quant_info[short_key]["rt stop"],
                        "y": formula_charge_to_quant_info[short_key]["amount"] / 2,
                        "text": "RT window stop",
                        "pos": 2,  # left,
                        "col": "blue",
                    },
                ]
            }

            for evidence_rt in formula_charge_to_quant_info[short_key]["evidence_rts"]:
                ablines[key].append({"v": evidence_rt, "lwd": 0.5, "col": "purple"})
                additional_legends[key].append(
                    {
                        "x": evidence_rt,
                        "y": 0,
                        "lwd": 0.5,
                        "col": "purple",
                        "text": "MS2 ident",
                        "pos": 4,
                        "srt": 45,  # rotate label
                    }
                )

            results.plot_MICs_2D(
                [key],
                file_name=None,
                rt_window=None,
                i_transform=None,
                xlimits=[
                    formula_charge_to_quant_info[short_key]["rt start"] - 0.05,
                    formula_charge_to_quant_info[short_key]["rt stop"] + 0.05,
                ],
                additional_legends=additional_legends,
                title=None,
                zlimits=None,
                ablines=ablines,
                graphics=graphics,
            )
            print("Plottted {0}".format(file_name))

    return


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
    else:
        main(ident_file=sys.argv[1], mzml_file=sys.argv[2])
