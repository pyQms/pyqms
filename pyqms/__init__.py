#!/usr/bin/env python
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
import os
from pkg_resources import parse_version

# version_info  = (0, 5, 0, 'beta')
# version = '0.5.0-beta'

__version_str__ = (
    open(os.path.join(os.path.dirname(__file__), "version.txt")).readline().strip()
)
__version__ = parse_version(__version_str__)

if not hasattr(sys, "version_info") or sys.version_info < (3, 5):
    raise RuntimeError("pyQms requires Python 3.5 or later.")

from . import knowledge_base
from .isotopologue_library import IsotopologueLibrary
from .results import Results
from .results import match
from .params import params

del sys
