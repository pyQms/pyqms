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
from urllib import request as request
import os
import shutil


def main():
    """
    Downloads the BSA.mzML example file also used in openMS and Ursgal.

    This file is ideally suited for the use with the example scripts to test
    pyQms.

    Will download the BSA1.mzML to the data folder ( ../data/BSA1.mzML )

    Usage:

        ./get_example_BSA_file.py

    """
    mzML_file = os.path.join(os.pardir, "data", "BSA1.mzML")
    if os.path.exists(mzML_file) is False:
        http_url = "http://sourceforge.net/p/open-ms/code/HEAD/tree/OpenMS/share/OpenMS/examples/BSA/BSA1.mzML?format=raw"
        basename = os.path.basename(http_url).replace("?", "")  # Win compatible
        output_path = os.path.join(os.path.dirname(mzML_file), basename)
        with open(output_path, "wb") as ooo:
            local_filename, headers = request.urlretrieve(
                http_url, filename=output_path
            )
        try:
            shutil.move("{0}?format=raw".format(mzML_file), mzML_file)
        except:
            shutil.move("{0}format=raw".format(mzML_file), mzML_file)
        print("Saved file as {0}".format(mzML_file))

    return


if __name__ == "__main__":
    main()
