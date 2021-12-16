#!/usr/bin/env python3.4
# encoding: utf-8

import glob
import os

# import pyqms.params
import pyqms
from pyqms.params import params_descriptions as params_descriptions
import pprint

if __name__ == "__main__":
    print(
        """
        Formatting params into rst files for the docs
"""
    )

    with open("source/params.rst", "w") as o:
        print(
            """
.. params:

.. default-domain:: py

.. _parameter section:

##########
Parameters
##########

pyQms default params, parsed from current params.py file.

.. note:: This sphinx source file was **auto-generated** using
    pyqms/docs/parse_params_for_docu.py, which parses pyqms/params.py
    Please **do not** modify this file directly, but rather the original
    parameter files!


""",
            file=o,
        )

        print(""".. code-block:: python\n""", file=o)
        print(">>> params = {", file=o)
        for pos, (k, v) in enumerate(sorted(pyqms.params.items())):
            print("\t'{0}'' : {1},".format(k, v), file=o)
        print("}", file=o)

        print(
            """

Descriptions
============

""",
            file=o,
        )
        for params_class in params_descriptions.keys():
            #             print("""
            # {0}
            # {1}
            #                 """.format(params_class,'-' * (len(params_class) + 1)),file=o)
            for v in params_descriptions[params_class]:
                print(
                    """
{0}
{1}

{2}

Default value: {3}
                    """.format(
                        v["key"],
                        '"' * (len(v["key"]) + 1),
                        v["description"],
                        v["default"],
                    ),
                    file=o,
                )
