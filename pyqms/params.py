#! /usr/bin/env python
# -*- coding: utf-8 -*-
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
# """
# The params holds all parameters required for pyqms.

# """

from collections import OrderedDict as odict

params = {
    "PERCENTILE_FORMAT_STRING": "{0:.3f}",
    "M_SCORE_THRESHOLD": 0.5,
    "ELEMENT_MIN_ABUNDANCE": 1e-3,
    "MIN_REL_PEAK_INTENSITY_FOR_MATCHING": 0.01,
    # which relative intensity is considered during matching ...
    "REQUIRED_PERCENTILE_PEAK_OVERLAP": 0.50,
    "MINIMUM_NUMBER_OF_MATCHED_ISOTOPOLOGUES": 2,
    "INTENSITY_TRANSFORMATION_FACTOR": 1e5,
    "UPPER_MZ_LIMIT": 2000,
    "LOWER_MZ_LIMIT": 150,
    "MZ_TRANSFORMATION_FACTOR": 10000,
    # ^- this value will be multiplied with mz values before converted
    # to integer. This means that values with a difference of
    # 0.1 ppm @ 1000 m/z won't be distinguishable
    "REL_MZ_RANGE": 5e-6,
    "REL_I_RANGE": 0.2,
    # ^- this is the rel intensity error range at rel. abundance of 1
    # lower abundance are scaled that their rel. abundance + i_range
    # sum up to the same value.
    # The lower the intensity the higher the higher the
    # imprecision ... we might just include it into the code
    # and leave out the params option ...
    "INTERNAL_PRECISION": 1000,
    "MAX_MOLECULES_PER_MATCH_BIN": 20,
    "MZ_SCORE_PERCENTILE": 0.4,
    # Intensity Score complements
    "SILAC_AAS_LOCKED_IN_EXPERIMENT": None,
    "BUILD_RESULT_INDEX": True,
    "MACHINE_OFFSET_IN_PPM": 0.0,
    # ^-- this will only be applied on calculated on mz values! not mass !
    "FIXED_LABEL_ISOTOPE_ENRICHMENT_LEVELS": {"15N": 0.994, "13C": 0.996, "2H": 0.994},
    "COLORS": {
        0.0: (37, 37, 37),
        0.1: (99, 99, 99),
        0.2: (150, 150, 150),
        0.3: (204, 204, 204),
        0.4: (247, 247, 247),
        0.5: (203, 27, 29),
        0.6: (248, 120, 72),
        0.7: (253, 219, 121),
        0.8: (209, 239, 121),
        0.9: (129, 202, 78),
        1: (27, 137, 62),
    },
}

params_descriptions = odict(
    [
        (
            "matching_and_scoring",
            [
                {
                    # 'simple_name' : 'Required percentile peak overlap between matched isotpologues and calculated isotopologues above threshold',
                    "description": """Defines the percentile how many theoretical
and measured peaks must overlap so that the match is considered further.
E.g. 0.5 dictates, that 2 of 4 peaks must ovelap""",
                    "default": params["REQUIRED_PERCENTILE_PEAK_OVERLAP"],
                    "key": "REQUIRED_PERCENTILE_PEAK_OVERLAP",
                },
                {
                    # 'simple_name' : 'Element minimum abundance',
                    "description": """Defines the minimum abundance of an element
to be considered for the calculation of the isotopologue(s)""",
                    "default": params["ELEMENT_MIN_ABUNDANCE"],
                    "key": "ELEMENT_MIN_ABUNDANCE",
                },
                {
                    # 'simple_name' : 'Min relative peak intensity required for matching',
                    "description": """Defines the relative minimum peak intensity
within an isotopologue to be considered for matching""",
                    "default": params["MIN_REL_PEAK_INTENSITY_FOR_MATCHING"],
                    "key": "MIN_REL_PEAK_INTENSITY_FOR_MATCHING",
                },
                {
                    # 'simple_name' : 'Relative intensity range',
                    "description": """Defines the relative intensity error range.
Represents the relative error to the most intense peak.""",
                    "default": params["REL_I_RANGE"],
                    "key": "REL_I_RANGE",
                },
                {
                    # 'simple_name' : 'Relative m/z range',
                    "description": """Defines the relative m/z error range or the
measuring precision of the used mass spectrometer. Is equal to the precision of
the used machine in parts per million (ppm)""",
                    "default": params["REL_MZ_RANGE"],
                    "key": "REL_MZ_RANGE",
                },
                {
                    # 'simple_name' : 'm/z score percentile',
                    "description": """Defines the weighting between the m/z error
and the intensity error for the total score. This weighting can be adjusted for
different mass spectrometers, depending on whether m/z or intensity can be
measured more accurately""",
                    "default": params["MZ_SCORE_PERCENTILE"],
                    "key": "MZ_SCORE_PERCENTILE",
                },
                {
                    # 'simple_name': 'Minimum number of isotopolgue matches required',
                    "description": """Number of isotopologue peaks that are required
to yield a mScore. Very small molecules may yield only one isotope peak
(monoisotopic peak) or the non-monoisotopic peaks have a very low abundance, so
that they ware not considered for macthing""",
                    "default": params["MINIMUM_NUMBER_OF_MATCHED_ISOTOPOLOGUES"],
                    "key": "MINIMUM_NUMBER_OF_MATCHED_ISOTOPOLOGUES",
                },
            ],
        ),
        (
            "measurement_and_reporting",
            [
                {
                    # 'simple_name' : 'Upper m/z limit',
                    "description": """Defines the maximum m/z value to be
considered by pyQms. Can be adjusted for better performance of pyQms or to
limit for the measuring range of the used mass spectrometer""",
                    "default": params["UPPER_MZ_LIMIT"],
                    "key": "UPPER_MZ_LIMIT",
                },
                {
                    # 'simple_name' : 'Lower m/z limit',
                    "description": """Defines the minimum m/z value to be
considered by pyQms. Can be adjusted for better performance of pyQms or to
limit for the measuring range of the used mass spectrometer""",
                    "default": params["LOWER_MZ_LIMIT"],
                    "key": "LOWER_MZ_LIMIT",
                },
                {
                    # 'simple_name' : 'Machine offset in ppm',
                    "description": """A mass spectrometer measuring error (constant
machine/calibration dependent mass or m/z offset) can be defined here in parts
per million (ppm)""",
                    "default": params["MACHINE_OFFSET_IN_PPM"],
                    "key": "MACHINE_OFFSET_IN_PPM",
                },
                {
                    # 'simple_name' : 'mScore threshold',
                    "description": """The minimum mScore, which should be reported.
Typically a mScore above 0.7 yields a FDR below 1%. Lower mScore thresholds
can be used to check for machine errors or to optimize matching of pulse-chase
samples""",
                    "default": params["M_SCORE_THRESHOLD"],
                    "key": "M_SCORE_THRESHOLD",
                },
                {
                    # 'simple_name' : 'Silac amino acids locked in experiment',
                    "description": """These aminoacids have always the defined
fixed SILCA modification and their atoms are not considered when calculating a
partially labeling percentile""",
                    "default": params["SILAC_AAS_LOCKED_IN_EXPERIMENT"],
                    "key": "SILAC_AAS_LOCKED_IN_EXPERIMENT",
                },
            ],
        ),
        (
            "internal",
            [
                {
                    # 'simple_name' : 'Percentile format string',
                    "description": """Defines the standard format string when
formatting labeling percentile float. Standard format considers three floating
points""",
                    "default": params["PERCENTILE_FORMAT_STRING"],
                    "key": "PERCENTILE_FORMAT_STRING",
                },
                {
                    # 'simple_name' : 'Internal precision',
                    "description": """Defines the internal precision for float to
int conversion""",
                    "default": params["INTERNAL_PRECISION"],
                    "key": "INTERNAL_PRECISION",
                },
                {
                    # 'simple_name' : 'Maximum molecules per match bin',
                    "description": """Defines the number of molecules per match bin.
Influences the matching speed""",
                    "default": params["MAX_MOLECULES_PER_MATCH_BIN"],
                    "key": "MAX_MOLECULES_PER_MATCH_BIN",
                },
                {
                    # 'simple_name' : 'm/z transformation factor',
                    "description": """All m/z values are transformed by this factor
This value will be multiplied with m/z values before converted to integer. This
means that values with a difference of 0.1 ppm @ 1000 m/z won't be
distinguishable""",
                    "default": params["MZ_TRANSFORMATION_FACTOR"],
                    "key": "MZ_TRANSFORMATION_FACTOR",
                },
                {
                    # 'simple_name' : 'Intensity transformation factor',
                    "description": """All intensities are transformed with this
factor""",
                    "default": params["INTENSITY_TRANSFORMATION_FACTOR"],
                    "key": "INTENSITY_TRANSFORMATION_FACTOR",
                },
                {
                    # 'simple_name' : 'Build result index',
                    "description": """The results are indexed for faster access""",
                    "default": params["BUILD_RESULT_INDEX"],
                    "key": "BUILD_RESULT_INDEX",
                },
            ],
        ),
    ]
)
