#!/usr/bin/env python
# encoding: utf-8

import pyqms
import sys
from pyqms.adaptors import calc_amount_function as calc_amount_function

TESTS = [
    {"input": {"rt": [], "i": [], "scores": [], "spec_ids": []}, "output": None},
    {
        "input": {
            "rt": [1, 2, 3],
            "i": [0, 10, 100],
            "scores": [0.7, 0.8, 0.9],
            "spec_ids": [1, 2, 3],
        },
        "output": {
            "max I in window": 100,
            "max I in window (rt)": 3,
            "max I in window (score)": 0.9,
            "sum I in window": 110,
            "auc in window": 60,
        },
    },
    {
        "input": {
            "rt": [1, 2, 3, 4],
            "i": [10, 100, 100, 10],
            "scores": [0.7, 0.8, 0.9, 0.8],
            "spec_ids": [1, 2, 3, 4],
        },
        "output": {
            "max I in window": 100,
            "max I in window (rt)": 2,
            "max I in window (score)": 0.8,
            "sum I in window": 220,
            "auc in window": 210,
        },
    },
    {
        "input": {
            "rt": [1, 2, 3, 4],
            "i": [10, 100, 120, 10],
            "scores": [0.7, 0.8, 0.9, 0.8],
            "spec_ids": [1, 2, 3, 4],
        },
        "output": {
            "max I in window": 120,
            "max I in window (rt)": 3,
            "max I in window (score)": 0.9,
            "sum I in window": 240,
            "auc in window": 230,
        },
    },
    {
        "input": {
            "rt": [1, 2, 3, 4],
            "i": [10, 10, 10, 10],
            "scores": [0.8, 0.8, 0.8, 0.8],
            "spec_ids": [1, 2, 3, 4],
        },
        "output": {
            "max I in window": 10,
            "max I in window (rt)": 1,
            "max I in window (score)": 0.8,
            "sum I in window": 40,
            "auc in window": 30,
        },
    },
]


def amount_test():
    for test_dict in TESTS:
        yield check_amount, test_dict


def check_amount(test_dict):

    result_dict = calc_amount_function(test_dict["input"])
    print(result_dict)
    if result_dict is None:
        assert test_dict["output"] is result_dict
    else:
        for key, value in result_dict.items():
            assert test_dict["output"][key] == value


if __name__ == "__main__":
    for test_dict in TESTS:
        check_amount(test_dict)
