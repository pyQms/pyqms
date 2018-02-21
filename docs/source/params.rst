
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



.. code-block:: python

>>> params = {
	'BUILD_RESULT_INDEX'' : True,
	'COLORS'' : {0.0: (37, 37, 37), 0.1: (99, 99, 99), 0.3: (204, 204, 204), 0.6: (248, 120, 72), 0.8: (209, 239, 121), 0.5: (203, 27, 29), 1: (27, 137, 62), 0.7: (253, 219, 121), 0.2: (150, 150, 150), 0.9: (129, 202, 78), 0.4: (247, 247, 247)},
	'ELEMENT_MIN_ABUNDANCE'' : 0.001,
	'FIXED_LABEL_ISOTOPE_ENRICHMENT_LEVELS'' : {'2H': 0.994, '13C': 0.996, '15N': 0.994},
	'INTENSITY_TRANSFORMATION_FACTOR'' : 100000.0,
	'INTERNAL_PRECISION'' : 1000,
	'LOWER_MZ_LIMIT'' : 150,
	'MACHINE_OFFSET_IN_PPM'' : 0.0,
	'MAX_MOLECULES_PER_MATCH_BIN'' : 20,
	'MINIMUM_NUMBER_OF_MATCHED_ISOTOPOLOGUES'' : 2,
	'MIN_REL_PEAK_INTENSITY_FOR_MATCHING'' : 0.01,
	'MZ_SCORE_PERCENTILE'' : 0.4,
	'MZ_TRANSFORMATION_FACTOR'' : 10000,
	'M_SCORE_THRESHOLD'' : 0.5,
	'PERCENTILE_FORMAT_STRING'' : {0:.3f},
	'REL_I_RANGE'' : 0.2,
	'REL_MZ_RANGE'' : 5e-06,
	'REQUIRED_PERCENTILE_PEAK_OVERLAP'' : 0.5,
	'SILAC_AAS_LOCKED_IN_EXPERIMENT'' : None,
	'UPPER_MZ_LIMIT'' : 2000,
}


Descriptions
============



REQUIRED_PERCENTILE_PEAK_OVERLAP
"""""""""""""""""""""""""""""""""

Defines the percentile how many theoretical
and measured peaks must overlap so that the match is considered further.
E.g. 0.5 dictates, that 2 of 4 peaks must ovelap

Default value: 0.5
                    

ELEMENT_MIN_ABUNDANCE
""""""""""""""""""""""

Defines the minimum abundance of an element
to be considered for the calculation of the isotopologue(s)

Default value: 0.001
                    

MIN_REL_PEAK_INTENSITY_FOR_MATCHING
""""""""""""""""""""""""""""""""""""

Defines the relative minimum peak intensity
within an isotopologue to be considered for matching

Default value: 0.01
                    

REL_I_RANGE
""""""""""""

Defines the relative intensity error range.
Represents the relative error to the most intense peak.

Default value: 0.2
                    

REL_MZ_RANGE
"""""""""""""

Defines the relative m/z error range or the
measuring precision of the used mass spectrometer. Is equal to the precision of
the used machine in parts per million (ppm)

Default value: 5e-06
                    

MZ_SCORE_PERCENTILE
""""""""""""""""""""

Defines the weighting between the m/z error
and the intensity error for the total score. This weighting can be adjusted for
different mass spectrometers, depending on whether m/z or intensity can be
measured more accurately

Default value: 0.4
                    

MINIMUM_NUMBER_OF_MATCHED_ISOTOPOLOGUES
""""""""""""""""""""""""""""""""""""""""

Number of isotopologue peaks that are required
to yield a mScore. Very small molecules may yield only one isotope peak
(monoisotopic peak) or the non-monoisotopic peaks have a very low abundance, so
that they ware not considered for macthing

Default value: 2
                    

UPPER_MZ_LIMIT
"""""""""""""""

Defines the maximum m/z value to be
considered by pyQms. Can be adjusted for better performance of pyQms or to
limit for the measuring range of the used mass spectrometer

Default value: 2000
                    

LOWER_MZ_LIMIT
"""""""""""""""

Defines the minimum m/z value to be
considered by pyQms. Can be adjusted for better performance of pyQms or to
limit for the measuring range of the used mass spectrometer

Default value: 150
                    

MACHINE_OFFSET_IN_PPM
""""""""""""""""""""""

A mass spectrometer measuring error (constant
machine/calibration dependent mass or m/z offset) can be defined here in parts
per million (ppm)

Default value: 0.0
                    

M_SCORE_THRESHOLD
""""""""""""""""""

The minimum mScore, which should be reported.
Typically a mScore above 0.7 yields a FDR below 1%. Lower mScore thresholds
can be used to check for machine errors or to optimize matching of pulse-chase
samples

Default value: 0.5
                    

SILAC_AAS_LOCKED_IN_EXPERIMENT
"""""""""""""""""""""""""""""""

These aminoacids have always the defined
fixed SILCA modification and their atoms are not considered when calculating a
partially labeling percentile

Default value: None
                    

PERCENTILE_FORMAT_STRING
"""""""""""""""""""""""""

Defines the standard format string when
formatting labeling percentile float. Standard format considers three floating
points

Default value: {0:.3f}
                    

INTERNAL_PRECISION
"""""""""""""""""""

Defines the internal precision for float to
int conversion

Default value: 1000
                    

MAX_MOLECULES_PER_MATCH_BIN
""""""""""""""""""""""""""""

Defines the number of molecules per match bin.
Influences the matching speed

Default value: 20
                    

MZ_TRANSFORMATION_FACTOR
"""""""""""""""""""""""""

All m/z values are transformed by this factor
This value will be multiplied with m/z values before converted to integer. This
means that values with a difference of 0.1 ppm @ 1000 m/z won't be
distinguishable

Default value: 10000
                    

INTENSITY_TRANSFORMATION_FACTOR
""""""""""""""""""""""""""""""""

All intensities are transformed with this
factor

Default value: 100000.0
                    

BUILD_RESULT_INDEX
"""""""""""""""""""

The results are indexed for faster access

Default value: True
                    
