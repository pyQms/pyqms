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
import sys
import pyqms


def main(args):
    """
    Uses a given peptide and charge and calculates and outputs the isotope
    envelope. Further options include metabolic labels and fixed labels.

    usage:

        ./view_isotopologue_overview.py molecule [charge [metabolic labels] [fixed labels]]

    e.g.:

        ./view_isotopologue_overview.py EILCEWRRAR 3 "{'15N' :[0,0.1]}" "{'R' :['C(-6) 13C(6)',''],'C':['C(1)O(2)','']}"

    Minimally a peptide and charge is required!

    """

    molecule = sys.argv[1]
    charges = [1]
    metabolic_labels = None
    fixed_labels = None
    if len(sys.argv) >= 3:
        charges = [int(sys.argv[2])]
        if len(sys.argv) >= 4:
            metabolic_labels = eval(sys.argv[3])
            if len(sys.argv) >= 5:
                fixed_labels = eval(sys.argv[4])

    lib = pyqms.IsotopologueLibrary(
        molecules=[molecule],
        charges=charges,
        metabolic_labels=metabolic_labels,
        fixed_labels=fixed_labels,
        params={"LOWER_MZ_LIMIT": 0},
    )

    for formula in lib.keys():
        for charge in charges:
            lib.print_overview(formula, charge=charge)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
    else:
        main(sys.argv)
