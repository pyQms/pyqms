#!/usr/bin/env python
# encoding: utf-8
"""
"""
import sys
import pyqms

MATCHED_PEAKS = {
    # mmz, mi, ri, cmz, ci in matched_peaks
    # expected score, expected scaling factor
    (1.0, 1.0): [(1, 1, 0.9, 1, 1), (0.5, 1, 0.2, 0.5, 1)],
    # scaling and score perfect
    (1.0, 2.0): [(1, 2, 0.9, 1, 1), (0.5, 2, 0.2, 0.5, 1)],
    # scaling times 2
    (1.0, 3.0): [(1, 3, 0.9, 1, 1), (0.123, 123, 0.0, 0.5, 1)],
    # scaling times 3 and last peak is ignored with scaling of 0.
    (0.0, 0.0): [(1.1, 3, 0.9, 1, 1), (0.55, 3, 0.2, 0.5, 1)],
    # bad match
    (0.5, 0.0): [(1, 0, 0.9, 1, 1), (0.5, 0, 0.2, 0.5, 1)],
    # perfect mz match and no i match, thus no scaling and no i score
    (0.5, 5.0): [(1, 5, 0.9, 0.5, 1), (0.5, 5, 0.2, 0.2, 1)],
    # # perfect i match and no mz match, thus scaling but no mz score
    (1.0, 6.0): [(1, 6, 1.0, 1, 1), (0.5, 6, 1.0, 0.5, 1), (None, None, 0.0, 0.5, 0.2)],
    # scaling and score perfect - None are ignored
}


LIB = pyqms.IsotopologueLibrary(
    molecules=["ELVISLIVES"],
    charges=[2],
    verbose=False,
    params={"MZ_SCORE_PERCENTILE": 0.5},
)


def score_test():
    for (expected_score, expected_scaling), matched_peaks in MATCHED_PEAKS.items():
        yield check_score, expected_score, matched_peaks


def scaling_test():
    for (expected_score, expected_scaling), matched_peaks in MATCHED_PEAKS.items():
        yield check_scaling, expected_scaling, matched_peaks


def check_scaling(expected_scaling, matched_peaks):
    score, scaling_factor = LIB.score_matches(matched_peaks, 0.5)
    assert float(expected_scaling) - scaling_factor <= sys.float_info.epsilon


def check_score(expected_score, matched_peaks):
    score, scaling_factor = LIB.score_matches(matched_peaks, 0.5)
    assert float(expected_score) - score <= sys.float_info.epsilon


if __name__ == "__main__":
    for (expected_score, expected_scaling), matched_peaks in MATCHED_PEAKS.items():
        # check_scaling(expected_scaling, matched_peaks)
        # check_score(expected_score, matched_peaks)
        score, scaling_factor = LIB.score_matches(matched_peaks, 0.5)
        print(score, scaling_factor, "Expected:", expected_score, expected_scaling)
