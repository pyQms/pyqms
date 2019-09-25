#!/usr/bin/env python
# encoding: utf-8

import pyqms

R = pyqms.Results()

VARIANTS = {
    # ( (1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2), 5 ) : \
    #     (1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5)
    # ,
    # ( (1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2), 1 ) : \
    #     (1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2)
    # ,
    ((0, 0, 0, 0, 0, 0, 0), 1): (0, 0, 0, 0, 0, 0, 0),
    ((0, 1, 2, 3, 2, 1, 0), 1): (0, 1, 2, 3, 2, 1, 0),
    ((0, 1, 2, 3, 2, 1, 0), 2): (0.5, 1.5, 3.0, 3.5, 3.0, 1.5, 0.5),
    ((2, 1, 2, 1, 2, 1, 2), 3): (
        3 / 3.0,
        5 / 3.0,
        4 / 3.0,
        5 / 3.0,
        4 / 3.0,
        5 / 3.0,
        3 / 3.0,
    ),
    ((0, 1, 2, 3, 2, 1, 0), 5): (
        3 / 5.0,
        6 / 5.0,
        8 / 5.0,
        9 / 5.0,
        8 / 5.0,
        6 / 5.0,
        3 / 5.0,
    ),
}


def smooth_test():
    for (input_list, k), expected_output in VARIANTS.items():
        yield smooth_list, input_list, k, expected_output


def smooth_list(input_list, k, expected_output):
    smoothed_list = R._smooth_list(input_list, k=k)
    # print('k', k)
    # print('input', input_list)
    # print('smoothend', tuple(smoothed_list))
    # print('Expected',expected_output)
    assert tuple(smoothed_list) == expected_output
    # assert pairs == expected_pairs


if __name__ == "__main__":
    for (input_list, k), expected_output in VARIANTS.items():
        smooth_list(input_list, k, expected_output)
