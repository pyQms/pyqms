#!/usr/bin/env python
# encoding: utf-8
import pyqms
import unittest
import math
import os
import pickle
import csv
from pyqms.adaptors import read_xlsx_file as read_xlsx_file
import pandas as pd


class TestResults(unittest.TestCase):
    def setUp(self):
        self.results = pyqms.Results(
            lookup={
                "formula to molecule": {
                    "C(37)H(59)N(9)O(16)": ["DDSPDLPK"],
                    "C(43)H(75)N(15)O(17)S(2)": [
                        "CCTESLVNR#Carbamidomethyl:1;Carbamidomethyl:2"
                    ],
                },
                "molecule to formula": {
                    "DDSPDLPK": "C(37)H(59)N(9)O(16)",
                    "CCTESLVNR#Carbamidomethyl:1;Carbamidomethyl:2": "C(43)H(75)N(15)O(17)S(2)",
                },
            }
        )
        keys = [
            ("BSA1.mzML", "C(37)H(59)N(9)O(16)", 2, (("N", "0.000"),)),
            ("BSA1.mzML", "C(37)H(59)N(9)O(16)", 2, (("N", "0.000"),)),
            ("BSA2.mzML", "C(43)H(75)N(15)O(17)S(2)", 3, (("N", "0.010"),)),
        ]
        values = [
            # mmz,         mi, rel_icmz,        ci
            (1337, 13.37, 1, 100, [(443.7112649, 100, 1, 443.7112649, 1)]),
            (1338, 13.38, 0.9, 100, [(443.7112649, 100, 1, 443.7112649, 1)]),
            (1337, 13.37, 1, 10, [(569.7526156, 10, 1, 569.7526156, 1)]),
        ]

        for key, value in zip(keys, values):
            self.results.add(key, value)

    def add_test(self):
        """
        Test the add funtion on a temporary result class
        Args:
            key (tuple):
                | file_name
                | formula
                | charge
                | label_percentiles
            value (tuple):
                | spec_id
                | rt
                | score
                | scaling_factor
                | peaks

        """
        tmp_results = pyqms.Results()
        keys = [("BSA1.mzML", "C(37)H(59)N(9)O(16)", 2, (("N", "0.000"),))]
        values = [(1337, 13.37, 1, 100, [(443.7112649, 100, 1, 443.7112649, 100)])]
        for key, value in zip(keys, values):
            tmp_results.add(key, value)
        assert len(tmp_results.index["files"]) == 1
        assert len(tmp_results.keys()) == 1
        return

    def parse_and_filter_test(self):
        """
        def _parse_and_filter(
            self,
            molecules         = None,
            charges           = None,
            file_names        = None,
            label_percentiles = None,
            formulas          = None
        ):
        """

        TESTS = [
            {
                "input": {"molecules": ["DDSPDLPK"]},
                "output": {"formula": "C(37)H(59)N(9)O(16)"},
            },
            {
                "input": {"charges": [3]},
                "output": {"formula": "C(43)H(75)N(15)O(17)S(2)"},
            },
            {
                "input": {"file_names": ["BSA1.mzML"]},
                "output": {"formula": "C(37)H(59)N(9)O(16)"},
            },
            {
                "input": {"formulas": ["C(37)H(59)N(9)O(16)"]},
                "output": {"formula": "C(37)H(59)N(9)O(16)"},
            },
            {
                "input": {"label_percentiles": [(("N", "0.010"),)]},
                "output": {"formula": "C(43)H(75)N(15)O(17)S(2)"},
            },
        ]
        for test_dict in TESTS:
            for key in self.results._parse_and_filter(**test_dict["input"]):
                assert key.formula == test_dict["output"]["formula"]

        return

    def extract_results_test(self):
        """
        extract_results(
            self,
            molecules         = None,
            charges           = None,
            file_names        = None,
            label_percentiles = None,
            formulas          = None,
            score_threshold   = None
        )

        keys    = [
            ('BSA1.mzML', 'C(37)H(59)N(9)O(16)', 2, (('N', '0.000'),) ),
            ('BSA2.mzML', 'C(43)H(75)N(15)O(17)S(2)', 3, (('N', '0.010'),) ),
        ]
        values  = [
            ( 1337, 13.37, 1, 100 , [(443.7112649, 100, 1, 443.7112649, 100 )] ),
            ( 1337, 13.37, 1, 10  , [(569.7526156, 10, 0.9, 569.7526156, 10 )] )
        ]

        """
        assert len(self.results.keys()) != 0
        TESTS = [
            {
                "input": {"molecules": ["DDSPDLPK"], "score_threshold": 0.95},
                "output": {
                    "formula": "C(37)H(59)N(9)O(16)",
                    "file_name": "BSA1.mzML",
                    "scaling_factor": 100,
                    "spec_id": 1337,
                },
            }
        ]
        for test_dict in TESTS:
            for key, n, entry in self.results.extract_results(**test_dict["input"]):
                print(key, entry)
                assert key.formula == test_dict["output"]["formula"]
                assert key.file_name == test_dict["output"]["file_name"]
                assert entry.scaling_factor == test_dict["output"]["scaling_factor"]
                assert entry.spec_id == test_dict["output"]["spec_id"]
            # print(self.results)
            # print(self.results.lookup)
            assert n == 0

    def extract_format_results_test(self):
        """
        format_results(
            self,
            molecules         = None,
            charges           = None,
            file_names        = None,
            label_percentiles = None,
            formulas          = None,
            score_threshold   = None
        )

        keys    = [
            ('BSA1.mzML', 'C(37)H(59)N(9)O(16)', 2, (('N', '0.000'),) ),
            ('BSA2.mzML', 'C(43)H(75)N(15)O(17)S(2)', 3, (('N', '0.010'),) ),
        ]
        values  = [
            ( 1337, 13.37, 1, 100 , [(443.7112649, 100, 1, 443.7112649, 100 )] ),
            ( 1338, 13.38, 0.9, 100 , [(443.7112649, 100, 1, 443.7112649, 1)] ),
            ( 1337, 13.37, 1, 10  , [(569.7526156, 10, 0.9, 569.7526156, 10 )] )
        ]

        """
        assert len(self.results.keys()) != 0
        TESTS = [
            {
                "output": [
                    {
                        "file_name": "BSA1.mzML",
                        "spec_id": 1337,
                        "formula": "C(37)H(59)N(9)O(16)",
                        "scaling_factor": 100,
                        "score": 1,
                        "charge": 2,
                    },
                    {
                        "file_name": "BSA1.mzML",
                        "spec_id": 1338,
                        "formula": "C(37)H(59)N(9)O(16)",
                        "scaling_factor": 100,
                        "score": 0.9,
                        "charge": 2,
                    },
                    {
                        "file_name": "BSA2.mzML",
                        "spec_id": 1337,
                        "formula": "C(43)H(75)N(15)O(17)S(2)",
                        "scaling_factor": 10,
                        "score": 1,
                        "charge": 3,
                    },
                ]
            }
        ]
        for test_dict in TESTS:
            values = self.results.format_all_results()

            assert isinstance(values, pd.DataFrame)

            for out_data in test_dict["output"]:
                result = values.loc[
                    (values["file_name"] == out_data["file_name"])
                    & (values["spec_id"] == out_data["spec_id"])
                ]
                assert (result["formula"] == out_data["formula"]).all()
                assert (result["scaling_factor"] == out_data["scaling_factor"]).all()
                assert (result["score"] == out_data["score"]).all()
                assert (result["charge"] == out_data["charge"]).all()

    def translate_test(self):
        """
        """
        assert self.results._translate_molecules_to_formulas(["DDSPDLPK"], None) == set(
            ["C(37)H(59)N(9)O(16)"]
        )

        assert self.results._translate_molecules_to_formulas(
            ["DDSPDLPK"], ["C(37)H(59)N(9)O(16)"]
        ) == set(["C(37)H(59)N(9)O(16)"])

        return

    def max_score_test(self):
        """
        max_score    = [0,None,None,None]
        max_score[0] = entry.score
        max_score[1] = key
        max_score[2] = i
        max_score[3] = entry
        """
        max_score_tuple = self.results.max_score(molecules=["DDSPDLPK"])
        assert max_score_tuple[0] == 1  # score
        assert max_score_tuple[3].scaling_factor == 100  # intensity

        assert self.results.max_score(molecules=["_DDSPDLPK_"]) == [0, None, None, None]
        return

    def max_i_amount_test(self):

        obj_for_calc_amount = {
            "rt": [1.0, 2.0, 3.0, 4.0],
            "i": [1, 100, 200, 1],
            "scores": [0.7, 0.8, 0.9, 1],
            "spec_ids": [1, 2, 3, 4],
        }

        max_I_dict = self.results.determine_max_itensity(obj_for_calc_amount)
        assert max_I_dict["max I in window"] == 200
        assert max_I_dict["max I in window (rt)"] == 3.0
        assert max_I_dict["max I in window (score)"] == 0.9
        return

    def measure_error_test(self):
        """
        We have perfect matches...
        """
        error_dict = self.results._determine_measured_error(
            score_threshold=0.5, plot=False
        )
        assert error_dict["mz_error"] == [0, 0, 0]
        assert error_dict["intensity_error"] == [0, 0, 0]

    def rpy2_import_test(self):
        """
        Test the import test function (test-ception)
        """
        try:
            import rpy2

            rpy2_present = True
        except:
            rpy2_present = False

        assert self.results._import_rpy2() is rpy2_present

        if rpy2_present:
            # R color part
            assert len(self.results._generate_r_colors("None", 10)) == 10

            rainbow = self.results._generate_r_colors("rainbow", 10)
            rainbow_r = self.results._generate_r_colors("rainbow_r", 10)
            rainbow_r.reverse()
            assert rainbow == rainbow_r

            plot_result_class = pickle.load(
                open(os.path.join("tests", "data", "test_BSA_pyqms_results.pkl"), "rb")
            )
            plot_name = os.path.join("tests", "data", "BSA_DDSPDLPK")
            for key in plot_result_class._parse_and_filter(molecules=["DDSPDLPK"]):
                # plot 3D
                plot_result_class.plot_MIC_3D(
                    key, file_name=plot_name, rt_window=None, i_transform=None
                )
                assert os.path.exists(plot_name + "_MIC_3D.pdf") is True

                # test fail
                plot_result_class.plot_MIC_3D(
                    key, file_name=plot_name, rt_window=[-2, -1], i_transform=None
                )
                # plot 2D
                graphics, grdevices = plot_result_class.init_r_plot(
                    plot_name + "_MIC_2D.pdf"
                )
                max_score_tuple = plot_result_class.max_score(molecules=["DDSPDLPK"])

                plot_result_class.plot_MICs_2D(
                    [key],
                    graphics=graphics,
                    rt_window=[28, 31],
                    ablines={
                        key: [
                            {
                                "v": max_score_tuple[3].rt,
                                "col": "gray",
                                "lty": "dashed",
                                "lwd": 0.4,
                            }
                        ]
                    },
                    additional_legends={
                        key: [{"text": "maxI RT", "x": max_score_tuple[3].rt, "y": 47}]
                    },
                )
                assert os.path.exists(plot_name + "_MIC_2D.pdf") is True

                # plot mz and i error function
                plot_result_class._determine_measured_error(
                    score_threshold=0.5,
                    filename=os.path.join(
                        "tests", "data", "mz_and_intensity_error_density_plot.pdf"
                    ),
                    plot=True,
                )

    def intensity_transformation_test(self):
        # normal way, no change
        i_label, i_transform_function = self.results._define_i_transformation()
        assert i_label == "Intensity [a.u.]"
        assert i_transform_function(10) == 10

        # log 2 transform, check
        i_label, i_transform_function = self.results._define_i_transformation(
            tag="log2"
        )
        assert i_label == "log2 Intensity [a.u.]"
        assert i_transform_function(10) == math.log(1 + 10, 2)

        # tag not defined, fall back to default
        i_label, i_transform_function = self.results._define_i_transformation(
            tag="log 2"
        )
        assert i_label == "Intensity [a.u.]"
        assert i_transform_function(10) == 10

        return

    def write_result_csv_test(self):

        output_file_name = os.path.join("tests", "data", "test_results.csv")
        self.results.lookup["formula to evidences"] = {
            "C(37)H(59)N(9)O(16)": {"DDSPDLPK": {"trivial_names": ["BSA"]}},
            "C(43)H(75)N(15)O(17)S(2)": {"CCTESLVNR": {"trivial_names": ["BSA"]}},
        }
        self.results.write_result_csv(output_file_name=output_file_name)
        assert os.path.exists(output_file_name) is True

        return

    def quant_summary_test(self):
        results_class = pickle.load(
            open(os.path.join("tests", "data", "test_BSA_pyqms_results.pkl"), "rb")
        )
        quant_summary_collector = []
        for pos, file_extension in enumerate(["xlsx", "csv", "extension_unknown"]):
            quant_summary_file = os.path.join(
                "tests", "data", "test_quant_summary.{0}".format(file_extension)
            )
            if pos == 2:
                with self.assertRaises(SystemExit) as system_exit_check:
                    results_class.calc_amounts_from_rt_info_file(
                        rt_info_file=quant_summary_file, rt_border_tolerance=None
                    )
                self.assertEqual(system_exit_check.exception.code, 1)
            else:
                # test csv writing
                results_class.write_rt_info_file(
                    output_file=quant_summary_file, rt_border_tolerance=1, update=True
                )
                assert os.path.exists(quant_summary_file) is True

                results_class.calc_amounts_from_rt_info_file(
                    rt_info_file=quant_summary_file,
                    rt_border_tolerance=1,
                    calc_amount_function=None,
                )
                if file_extension == "csv":
                    read_in_line_dicts = []
                    for line_dict in csv.DictReader(open(quant_summary_file, "r")):
                        read_in_line_dicts.append(line_dict)
                else:
                    read_in_line_dicts = read_xlsx_file(quant_summary_file)
                for line_dict in read_in_line_dicts:
                    assert line_dict["max I in window"] != ""
                    assert line_dict["max I in window (rt)"] != ""
                    assert line_dict["max I in window (score)"] != ""
                assert len(read_in_line_dicts) == len(results_class.keys())
                quant_summary_collector.append(read_in_line_dicts)
        assert len(quant_summary_collector[0]) == len(quant_summary_collector[1])
        return

    def label_percentile_test(self):
        list_of_xzy_dicts = self.results.determine_label_efficiency(element="N")
        assert len(list_of_xzy_dicts) == 3

        return

    def group_14N_15N_pairs_test(self):
        tmp_metabolic_labeling_class = pyqms.Results(
            params={"PERCENTILE_FORMAT_STRING": "{0:.3f}", "BUILD_RESULT_INDEX": True}
        )
        keys = [
            ("run1", "C(37)H(59)N(9)O(16)", 2, (("N", "0.000"),)),
            ("run1", "C(37)H(59)N(9)O(16)", 2, (("N", "0.990"),)),
        ]
        values = [
            (1337, 13.37, 1, 100, [(443.7112649, 100, 1, 443.7112649, 1)]),
            (1337, 13.37, 0.9, 200, [(452.7112649, 100, 1, 443.7112649, 1)]),
        ]
        for key, value in zip(keys, values):
            tmp_metabolic_labeling_class.add(key, value)
        all_key_pairs = []
        for light_key, heavy_key in tmp_metabolic_labeling_class._group_14N_15N_pairs():
            all_key_pairs.append((light_key, heavy_key))
        print(all_key_pairs)
        assert len(all_key_pairs) == 1
        # assert sorted(keys) == sorted(all_key_pairs)
        return

    def window_definition_and_curator_test(self):
        """
        {
            'RT'          : float(line_dict['Retention Time (s)']) / 60.0, # always in min
            'score'       : float(line_dict['PEP']),
            'score_field' : 'PEP',

        }

        """
        TESTS = [
            {
                "input": {
                    "evidence_dict": {
                        "peptide_1": {"evidences": [{"RT": 2}, {"RT": 4}]}
                    },
                    "rt_tolerance": 1,
                },
                "output": {"peptide_1": {"rt_window": [2, 4]}},
            },
            {
                "input": {
                    "evidence_dict": {
                        "peptide_1": {"evidences": [{"RT": 2}, {"RT": 4}]},
                        "peptide_2": {"evidences": [{"RT": 3}, {"RT": 5}]},
                    },
                    "rt_tolerance": 1,
                },
                "output": {
                    "peptide_1": {"rt_window": [2, 4], "upper_window_border": -0.5},
                    "peptide_2": {"rt_window": [3, 5], "lower_window_border": -0.5},
                },
            },
            {
                "input": {
                    "evidence_dict": {
                        "peptide_1": {"evidences": [{"RT": 2}, {"RT": 4}]},
                        "peptide_2": {"evidences": [{"RT": 3}, {"RT": 5}]},
                        "peptide_3": {"evidences": [{"RT": 4}, {"RT": 6}]},
                    },
                    "rt_tolerance": 1,
                },
                "output": {
                    "peptide_1": {"rt_window": [2, 4], "upper_window_border": -0.5},
                    "peptide_2": {
                        "rt_window": [3, 5],
                        "lower_window_border": -0.5,
                        "upper_window_border": -0.5,
                    },
                    "peptide_3": {"rt_window": [4, 6], "lower_window_border": -0.5},
                },
            },
            {
                "input": {
                    "evidence_dict": {
                        "peptide_1": {"evidences": [{"RT": 2}, {"RT": 3}]},
                        "peptide_2": {"evidences": [{"RT": 4}, {"RT": 5}]},
                    },
                    "rt_tolerance": 1,
                },
                "output": {
                    "peptide_1": {"rt_window": [2, 3], "upper_window_border": 0.5},
                    "peptide_2": {"rt_window": [4, 5], "lower_window_border": 0.5},
                },
            },
        ]

        keys_2_check = ["lower_window_border", "upper_window_border"]

        for test_dict in TESTS:
            output = self.results.curate_rt_windows(**test_dict["input"])
            for peptide in output.keys():
                assert (
                    output[peptide]["rt_window"]
                    == test_dict["output"][peptide]["rt_window"]
                )
                for key in keys_2_check:
                    if key in test_dict["output"][peptide].keys():
                        assert test_dict["output"][peptide][key] == output[peptide][key]
        return
