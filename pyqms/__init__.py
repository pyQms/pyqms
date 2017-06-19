#!/usr/bin/env python3.3
# encoding: utf-8
"""
    pyQms
    -----

    Python module for fast and accurate mass spectrometry data quantification

    :license: MIT, see LICENSE.txt for more details

    Authors:

        * Leufken, J.
        * Niehues, A.
        * Wessel, F.
        * Sarin, L.P.
        * Hippler, M.
        * Leidel, S.A.
        * Fufezan, C.

"""
from __future__ import absolute_import
import sys
version_info  = (0, 5, 0, 'beta')
version = '0.5.0-beta'

if not hasattr(sys, "version_info") or sys.version_info < (3, 3):
    raise RuntimeError("pyQms requires Python 3.3 or later.")

from . import knowledge_base
from .isotopologue_library import IsotopologueLibrary
from .chemical_composition import ChemicalComposition
from .unimod_mapper import UnimodMapper
from .results import Results
from .results import match
from .params import params
del sys
