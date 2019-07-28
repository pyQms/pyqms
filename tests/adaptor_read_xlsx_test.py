#!/usr/bin/env python
# encoding: utf-8
"""

adaptor test: read_xlsx_file

"""


from pyqms.adaptors import read_xlsx_file as read_xlsx_file
import os

TESTS = [
    {
        "input": os.path.join("tests", "data", "test_BSA_quant_summary.xlsx"),
        "output": {
            "molecules": ["DDSPDLPK", "CCTESLVNR#Carbamidomethyl:1;Carbamidomethyl:2"]
        },
    }
]


def adaptor_test():
    for test_dict in TESTS:
        yield adaptor_check, test_dict


def adaptor_check(test_dict):

    list_of_row_dicts = read_xlsx_file(test_dict["input"])

    results = {"molecules": list()}
    for row_dict in list_of_row_dicts:
        results["molecules"].append(row_dict["molecule"])

    assert len(results["molecules"]) == len(test_dict["output"]["molecules"])


if __name__ == "__main__":
    for test_dict in TESTS:
        adaptor_check(test_dict)
