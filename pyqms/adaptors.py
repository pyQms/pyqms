#! /usr/bin/env python
# -*- coding: utf-8 -*-
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
import os
from collections import defaultdict as ddict
import csv
import pyqms
from copy import deepcopy as dc
import sys
import codecs
import re

POST_EXPERIMENTAL_MODIFICATIONS = ["Carbamidomethyl"]

ELEMENT_REPLACEMENT_DICT = {"N": "14N"}
PARAM_TYPE_LOOKUP = {
    "PERCENTILE_FORMAT_STRING": str,
    "M_SCORE_THRESHOLD": float,
    "ELEMENT_MIN_ABUNDANCE": float,
    "MIN_REL_PEAK_INTENSITY_FOR_MATCHING": float,
    "REQUIRED_PERCENTILE_PEAK_OVERLAP": float,
    "MINIMUM_NUMBER_OF_MATCHED_ISOTOPOLOGUES": int,
    "INTENSITY_TRANSFORMATION_FACTOR": int,
    "UPPER_MZ_LIMIT": int,
    "LOWER_MZ_LIMIT": int,
    "MZ_TRANSFORMATION_FACTOR": int,
    "REL_MZ_RANGE": float,
    "REL_I_RANGE": float,
    "INTERNAL_PRECISION": int,
    "MAX_MOLECULES_PER_MATCH_BIN": int,
    "MZ_SCORE_PERCENTILE": float,
    "SILAC_AAS_LOCKED_IN_EXPERIMENT": str,
    "BUILD_RESULT_INDEX": bool,
    "MACHINE_OFFSET_IN_PPM": float,
}


def _parse_evidence_and_format_fixed_labels(data=None):
    """

    Reformats input params to pyQms compatible params. Additionally evidence
    files are read in and the fixed labels are reformatted (stripped
    from the modifications if peptides are read in) as required by pyQms.
    This is especially required if data/samples contains Carbamidomethylation
    as modification and the sample was e.g. 15N labeled. This ensures that the
    nitrogens pools of the peptides (which are 15N labeled) do not mix up with
    the nitrogen pool of the Carbamidomethylation (14N since intriduced during
    sample preparation).
    All fixed modifications needs to be specified so that is can be ignored from
    the input evidence file but correctly formatted for the parameters.


    Example format::

        {
            'molecules' : {'PEPTIDEA',...},
            'evidences' : {
                'C18H36O18N9' : {
                    'PEPTIDEA' : {
                        'evidences' : [
                            {
                                'RT':13.37,
                                'score':0.01,
                                'score_field':'PEP'
                            },
                            ...
                        ],
                        'trivial_names':[
                            'PROTEIN_NAME',
                            'PATHWAY_NAME',
                            ...
                        ]
                },
            },
            'charges'   : {1,2,...},
            'params'    : {
                'MACHINE_OFFSET_IN_PPM':0,
                ...
            },

        }


    Example of data passed::

        {
            'params': {'measurement_and_reporting': {'NAME': 'default'}},
            'fixed_labels': [
                {
                    'modification': {
                        'unimodID': '4',
                        'specificity_sites': ['C'],
                        'mono_mass': 57.021464,
                        'element': {'O': 1, 'H': 3, 'N': 1, 'C': 2},
                        'name': 'Carbamidomethyl'
                    },
                    'AA': 'C'
                }
            ],
            'molecules': 'AA',
            'metabolic_labels': [{'modification': '0, 0.99', 'atom': '15N'}],
            'charges': [1, 2, 3, 4, 5],
            'file': '/BSA1.mzML'
        }

        Returns:
        dict: molecules, evidences, correctly fomatted fixed labels, charges and
            parameters


    """
    r = {
        "molecules": set(),
        "evidences": None,
        "fixed_labels": None,
        "params": {},
        "charges": set(),
    }

    for charge in data["charges"]:

        if str(charge) in ["True", "False"]:
            continue
        try:
            int(charge)
            r["charges"].add(int(charge))
        except:
            pass

    for major_cat in data["params"].keys():
        for k, v in data["params"][major_cat].items():
            if k == "NAME":
                continue
            converted_value = PARAM_TYPE_LOOKUP[k](v)
            r["params"][k] = converted_value

    cc_factory = pyqms.chemical_composition.ChemicalComposition()
    tmp_fixed_labels = None
    if "fixed_labels" in data.keys() and len(data["fixed_labels"]) != 0:
        tmp_fixed_labels = {}
        for fl_dict in data["fixed_labels"]:
            if fl_dict["modification"]["name"] != "None":
                cc_factory.add_chemical_formula(fl_dict["modification"]["element"])

                if fl_dict["modification"]["name"] in POST_EXPERIMENTAL_MODIFICATIONS:
                    cc_factory["14N"] = cc_factory["N"]
                    del cc_factory["N"]

            aa = fl_dict["AA"]
            if aa not in tmp_fixed_labels.keys():
                tmp_fixed_labels[aa] = []
            tmp_fixed_labels[aa].append(
                {
                    "element_composition": dc(cc_factory),
                    "evidence_mod_name": fl_dict["modification"]["name"],
                }
            )
            cc_factory.clear()

    evidence_file = data.get("evidence_file", None)

    evidence_file_list = []
    if evidence_file is not None:
        evidence_file_list = [evidence_file]
    if "evidence_score_field" not in data.keys():
        data["evidence_score_field"] = "PEP"  #  default
    formatted_fixed_labels, evidence_lookup, molecule_list = parse_evidence(
        fixed_labels=tmp_fixed_labels,
        evidence_files=evidence_file_list,
        molecules=data["molecules"],
        evidence_score_field=data["evidence_score_field"],
    )
    r["fixed_labels"] = formatted_fixed_labels
    r["molecules"] = molecule_list
    r["evidences"] = evidence_lookup

    return r


def parse_evidence(
    fixed_labels=None,
    evidence_files=None,
    molecules=None,
    evidence_score_field=None,
    return_raw_csv_data=False,
):
    """
    Reads in the evidence file and returns the final formatted fixed labels,
    the evidence lookup, which is passed to the isotopologue library and the
    final formatted molecules (fixed labels are stripped form the molecules).

    Note:

        Output .csv files from `Ursgal`_ (`Documentation`_) can directly be
        used. Also `mzTab`_ files can be used as input.

    .. _Ursgal:
        https://github.com/ursgal/ursgal

    .. _Documentation:
        http://ursgal.readthedocs.io/en/latest/

    .. _mzTab:
        http://www.psidev.info/mztab

    Args:
        fixed_labels (dict): dict with fixed labels, example format is shown
            below.
        evidence_files (list): list of evidence file paths.
        molecules (list): list of additional molecules
        evidence_score_field (str): specify fieldname which holds the search
            engine score (Default is "PEP")

    Example fixed label format::

        {
            'C' : [
                {
                    'element': {
                        'O': 1,
                        'H': 3,
                        '14N': 1,
                        'C': 2
                    },
                    'evidence_mod_name': 'Carbamidomethyl'
                },
            ]
        }

    Returns:

        tuple: final formatted fixed label dict, evidence lookup, list of molecules

    """
    if molecules is None:
        molecules = []
    if evidence_score_field is None:
        evidence_score_field = "PEP"  #  default

    unimod_parser = pyqms.UnimodMapper()

    fixed_mod_lookup = {}
    amino_acid_2_fixed_mod_name = ddict(list)

    formatted_fixed_labels = None
    evidence_lookup = None
    molecule_set = set()

    all_fixed_mod_names = set()

    if fixed_labels is not None and len(fixed_labels.keys()) != 0:
        formatted_fixed_labels = {}
        for aa, fixed_mod_info_dict_list in fixed_labels.items():
            for fixed_mod_info_dict in fixed_mod_info_dict_list:
                if isinstance(fixed_mod_info_dict["element_composition"], dict):
                    tmp_cc_factory = pyqms.chemical_composition.ChemicalComposition()
                    tmp_cc_factory.add_chemical_formula(
                        fixed_mod_info_dict["element_composition"]
                    )
                else:
                    tmp_cc_factory = fixed_mod_info_dict["element_composition"]
                # print(type(tmp_cc_factory))
                # print(fixed_mod_info_dict)
                if aa not in formatted_fixed_labels.keys():
                    formatted_fixed_labels[aa] = []
                formatted_fixed_labels[aa].append(tmp_cc_factory.hill_notation_unimod())
                # save it under name and amino acid!
                fixed_mod_lookup[fixed_mod_info_dict["evidence_mod_name"]] = dc(
                    tmp_cc_factory
                )
                amino_acid_2_fixed_mod_name[aa].append(
                    fixed_mod_info_dict["evidence_mod_name"]
                )
                all_fixed_mod_names.add(fixed_mod_info_dict["evidence_mod_name"])
                tmp_cc_factory.clear()

    cc_factory = pyqms.chemical_composition.ChemicalComposition()

    # this is the lookup for the lib with the evidences
    # tmp_evidences = ddict(list)
    tmp_evidences = {}

    csv_raw_data_to_return = {}
    # tmp_charges_of_evidences = set()
    for evidence_file in evidence_files:
        input_is_csv = False
        evidence_lookup = {}
        with codecs.open(
            evidence_file, mode="r", encoding="utf-8"
        ) as openend_evidence_file:
            # first buffer the file here depending on mztab andf csv input
            if evidence_file.upper().endswith("CSV"):
                dict_reader = csv.DictReader(openend_evidence_file)
                modification_fieldname = "Modifications"
                rt_fieldname = "Retention Time (s)"
                seq_fieldname = "Sequence"
                input_is_csv = True
            elif evidence_file.upper().endswith("MZTAB"):
                dict_reader = csv.DictReader(
                    [row for row in openend_evidence_file if row[:3] in ["PSM", "PSH"]],
                    delimiter="\t",
                )
                modification_fieldname = "modifications"
                rt_fieldname = "retention_time"
                seq_fieldname = "sequence"
            else:
                print(
                    "The format {0} is not recognized by the pyQms adaptor function".format(
                        os.path.splitext(evidence_file)[1]
                    )
                )

            input_buffer = []
            for line_dict in dict_reader:
                input_buffer.append(line_dict)
            csv_raw_data_to_return[evidence_file] = input_buffer
            for line_dict in input_buffer:

                modifications = line_dict.get(modification_fieldname, "")
                if modifications == "":
                    molecule = line_dict[seq_fieldname]
                else:
                    if input_is_csv:
                        formatted_mods = line_dict[modification_fieldname]
                    else:
                        formatted_mods = []
                        # 2-UNIMOD:4,3-UNIMOD:4
                        for pos_and_unimod_id in line_dict[
                            modification_fieldname
                        ].split(","):
                            pos, unimod_id = pos_and_unimod_id.split("-")
                            unimod_name = unimod_parser.id2name(unimod_id.split(":")[1])
                            formatted_mods.append("{0}:{1}".format(unimod_name, pos))
                        formatted_mods = ";".join(formatted_mods)

                    molecule = "{0}#{1}".format(
                        line_dict[seq_fieldname], formatted_mods
                    )

                dict_2_append = {}
                rt = line_dict.get(rt_fieldname, "")
                # seconds is the standard also for mzTab
                if rt != "":
                    dict_2_append["RT"] = float(rt) / 60.0  # always in min

                score = line_dict.get(evidence_score_field, "")
                if score != "":
                    dict_2_append["score"] = float(score)
                    dict_2_append["score_field"] = evidence_score_field
                else:
                    dict_2_append["score"] = "None"
                    dict_2_append["score_field"] = "None"

                if molecule not in tmp_evidences.keys():
                    tmp_evidences[molecule] = {"evidences": [], "trivial_names": set()}
                for trivial_name_key in [
                    "proteinacc_start_stop_pre_post_;",  # old ursgal style
                    "trivial_name",  # self defined name
                    "Protein ID",  # new ursgal style
                    "accession",  # mzTab style
                ]:
                    additional_name = line_dict.get(trivial_name_key, "")
                    if additional_name != "":
                        # use set to remove double values
                        tmp_evidences[molecule]["trivial_names"].add(additional_name)
                        if 'trivial_name' not in dict_2_append.keys():
                            dict_2_append['trivial_name'] = additional_name
                        else:
                            dict_2_append['trivial_name'] += ';{0}'.format(additional_name)                
                tmp_evidences[molecule]["evidences"].append(dict_2_append)

    mod_pattern = re.compile(r""":(?P<pos>[0-9]*$)""")

    all_molecules = list(molecules)

    if len(tmp_evidences.keys()) > 0:
        all_molecules += list(tmp_evidences.keys())

    for molecule_and_mods in sorted(all_molecules):
        # try to convert trivial name set to list for conveniences
        try:
            tmp_evidences[molecule_and_mods]["trivial_names"] = sorted(
                list(set(tmp_evidences[molecule_and_mods]["trivial_names"]))
            )
        except:
            pass
        # print(molecule_and_mods)
        if "#" in molecule_and_mods:
            molecule, modifications = molecule_and_mods.split("#")
        else:
            molecule = molecule_and_mods
            modifications = None
        fixed_label_mod_addon_names = []
        if modifications is not None:
            mods_to_delete = []
            mod_list = modifications.split(";")
            for pos_in_mod_list, mod_and_pos in enumerate(mod_list):
                # OLD STYLE, no ':' in mod allowed!
                # mod, pos = mod_and_pos.split(':')
                # NEW STYLE, SILAC does not crash...
                for match in mod_pattern.finditer(mod_and_pos):
                    pos = int(match.group("pos"))
                    mod = mod_and_pos[: match.start()]
                    break

                modded_aa = molecule[int(pos) - 1]

                if (
                    formatted_fixed_labels is not None
                    and modded_aa in formatted_fixed_labels.keys()
                    and mod in all_fixed_mod_names
                ):
                    fixed_label_mod_addon_names.append(mod)
                    mods_to_delete.append(pos_in_mod_list)

            for modpos_2_remove in sorted(mods_to_delete, reverse=True):
                mod_list.pop(modpos_2_remove)

            if len(mod_list) > 0:
                molecule = "{0}#{1}".format(molecule, ";".join(mod_list))
            else:
                # nosetest does not line else and pass
                # molecule = molecule
                pass
        else:
            # fail check if fixed mod is not in the modifications!
            # add all fixed modification!
            if formatted_fixed_labels is not None:
                for aa in molecule:
                    if aa in formatted_fixed_labels.keys():
                        for mod_name in amino_acid_2_fixed_mod_name[aa]:
                            fixed_label_mod_addon_names.append(mod_name)
        # print(molecule)
        if molecule.startswith("+"):
            cc_factory.add_chemical_formula(molecule)
        else:
            cc_factory.use(molecule)
        if len(fixed_label_mod_addon_names) != 0:
            for fixed_mod_name in fixed_label_mod_addon_names:
                cc_factory.add_chemical_formula(fixed_mod_lookup[fixed_mod_name])
        complete_formula = cc_factory.hill_notation_unimod()

        molecule_set.add(molecule)
        if molecule_and_mods in tmp_evidences.keys():
            if complete_formula not in evidence_lookup.keys():
                evidence_lookup[complete_formula] = {}
            evidence_lookup[complete_formula][molecule_and_mods] = tmp_evidences[
                molecule_and_mods
            ]

        cc_factory.clear()

    molecule_list = list(molecule_set)

    if return_raw_csv_data:
        return (
            formatted_fixed_labels,
            evidence_lookup,
            molecule_list,
            csv_raw_data_to_return,
        )
    else:
        return formatted_fixed_labels, evidence_lookup, molecule_list


def calc_amount_function(obj_for_calc_amount):
    """
    Calculates actual molecule amounts. Three types of amounts are
    calculated for a matched isotope chromatogram (MIC), maximum intensity,
    summed up intensity and area under curve.
    Additionally the score and the retention time at the maximum
    intensity are determined.

    A test function exists to check correct amount determination.

    Returned keys in amound dict:

        * 'max I in window'
        * 'max I in window (rt)'
        * 'max I in window (score)'
        * 'auc in window'
        * 'sum I in window'

    Returns:

        dict: amount dict with keys shown above.

    """
    amount_dict = None
    if len(obj_for_calc_amount["i"]) != 0:
        amount_dict = {}
        maxI = max(obj_for_calc_amount["i"])
        index_of_maxI = obj_for_calc_amount["i"].index(maxI)
        amount_rt = obj_for_calc_amount["rt"][index_of_maxI]
        amount_score = obj_for_calc_amount["scores"][index_of_maxI]

        amount_dict["max I in window"] = maxI
        amount_dict["max I in window (rt)"] = amount_rt
        amount_dict["max I in window (score)"] = amount_score
        amount_dict["sum I in window"] = sum(obj_for_calc_amount["i"])
        amount_dict["auc in window"] = 0
        x0, y0, x1, y1 = 0, 0, 0, 0
        for pos_in_profile, intensity in enumerate(obj_for_calc_amount["i"]):
            if pos_in_profile == 0:
                continue
            # for auc calculation we need to be at position 2...
            x0 = obj_for_calc_amount["rt"][pos_in_profile - 1]
            # i.e. the last rt
            y0 = obj_for_calc_amount["i"][pos_in_profile - 1]
            # i.e. the last intensity
            x1 = obj_for_calc_amount["rt"][pos_in_profile]
            y1 = intensity

            xspace = x1 - x0
            height_of_triangle = abs(y1 - y0)
            square = xspace * y0
            triangle = 0.5 * (xspace * height_of_triangle)
            amount_dict["auc in window"] += square
            if y0 < y1:
                amount_dict["auc in window"] += triangle
            elif y0 > y1:
                amount_dict["auc in window"] -= triangle
            else:
                amount_dict["auc in window"] += 0

    return amount_dict


def read_xlsx_file(xlsx_file):
    list_of_row_dicts = []
    try:
        # from openpyxl import Workbook
        from openpyxl import load_workbook
    except:
        print("{0} is not installed, please install it and try again")
        print("pip3.4 install openpyxl")
        exit()
    wb = load_workbook(filename=xlsx_file)
    ws = wb.active
    tmp_file_headers = []
    for row_count, row in enumerate(ws.rows):
        if row_count == 0:
            for cell in row:
                tmp_file_headers.append(cell.value)
        else:
            row_dict = {}
            for cell_pos, cell in enumerate(row):
                if cell.value is None:
                    cell.value = ""
                row_dict[tmp_file_headers[cell_pos]] = cell.value
            list_of_row_dicts.append(row_dict)

    return list_of_row_dicts
