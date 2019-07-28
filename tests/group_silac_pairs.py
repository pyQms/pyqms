#!/usr/bin/env python
# encoding: utf-8

import pyqms
import sys

R = pyqms.Results()
R.lookup = {}
R.lookup["molecule fixed label variations"] = {
    "ORGINAL_MOLECULE": [
        "R0_C1_R0",  # you flip !
        "R0_C0_R0",
        "K0_C1_K0",  # you flip !
        "K1_C0_K1",  # you flip !
        "R0_C1_R1",  # decoys since mixed states
        "R1_C1_R1",
    ]
}

VARIANTS = {
    (("R0", "R1"),): [["R0_C0_R0", "R1_C0_R1"], ["R0_C1_R0", "R1_C1_R1"]],
    (("K0", "K1"),): [["K0_C1_K0", "K1_C1_K1"]],
    (("k0", "K1"),): [["K0_C1_K0", "K1_C1_K1"]],
    (("R0", "R1"), ("K0", "K1")): [
        ["R0_C0_R0", "R1_C0_R1"],
        ["R0_C1_R0", "R1_C1_R1"],
        ["K0_C1_K0", "K1_C1_K1"],
    ],
    (("R0", "R1"), ("K1", "K0")): [
        ["R0_C0_R0", "R1_C0_R1"],
        ["R0_C1_R0", "R1_C1_R1"],
        ["K1_C0_K1", "K0_C0_K0"],
    ],
}


def pair_silac_test():
    for silac_pairs, expected_pairs in VARIANTS.items():
        yield group_silac, silac_pairs, expected_pairs


def group_silac(silac_pairs, expected_pairs):
    grouped_silac_pairs = sorted(R._group_silac_pairs(silac_pairs=silac_pairs))
    for pos, pair in enumerate(sorted(expected_pairs)):
        print("_".join(pair), grouped_silac_pairs[pos])
        assert "_".join(pair) == "_".join(grouped_silac_pairs[pos])
    # assert pairs == expected_pairs


if __name__ == "__main__":
    for silac_pairs, expected_pairs in VARIANTS.items():
        group_silac(silac_pairs, expected_pairs)
