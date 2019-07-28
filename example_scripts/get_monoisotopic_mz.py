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
    Uses a given peptide and charge and returns the monoisotopic mz, i.e.
    postion 0 in the isotope envelope.

    usage:

        ./get_monoisotopic_mz.py <molecule> <charge>

    e.g.:

        ./get_monoisotopic_mz.py EILCEWRRAR 3

    """

    molecule = sys.argv[1]
    charge = int(sys.argv[2])

    lib = pyqms.IsotopologueLibrary(
        molecules=[molecule],
        charges=[charge],
        metabolic_labels=None,
        fixed_labels=None,
        verbose=False,
    )

    for formula in lib.keys():
        print(
            "Peptide {0} with formula {1} has a monoisotopic m/z of {2} @ charge {3}".format(
                molecule,
                formula,
                lib[formula]["env"][(("N", "0.000"),)][charge]["mz"][0],
                charge,
            )
        )
    return


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
    else:
        main(sys.argv)
