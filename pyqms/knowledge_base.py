#!/usr/bin/env python3.4
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
# The knowledge base holds all constants required for pyqms.

# ISOTOPIC_DISTRIBUTIONS correspond to Berglund M. & Wieser M.E.
# Isotopic compositions of the elements 2009 (IUPAC Technical Report)
# Pure Appl. Chem., 2011, Vol. 83, No. 2, pp. 397-410
# http://dx.doi.org/10.1351/PAC-REP-10-06-02
# Published online 2011-01-14
# http://www.ciaaw.org/pubs/TICE2009.pdf


# """

isotopic_distributions = {
    #               .. will be sorted while importing ..
    # Element:  [ (mass1, relAbundance1), ..., (massN, relAbundanceN)]
    #
    # NOTE: This variable will be update if fixed_labels is set.
    #       Ultimately each Isotopoloque class instant will store
    #       the modified dictionary under self.isotopic_distributions.
    #
    # NOTE: lib._extend_kb_with_fixed_labels treats list[0] element
    #       as natural isotope. This is not true for all elements
    #
    # Isotopic masses verified using
    # http://ciaaw.org/atomic-masses.htm
    # Wang, M., Audi, G., Wapstra, A. H., Kondev, F. G., MacCormick, M.,
    # Xu, X. and Pfeiffer, B. (2012)
    # The Ame2012 atomic mass evaluation. Chinese Phys. C 36, 1603–2014
    'H': [ ( 1.0078250322, 0.999885), ( 2.0141017781, 0.000115 ) ],
    'C': [ (12.0000,       0.9893  ), (13.003354835 , 0.0107   ) ],
    'N': [ (14.003074004,  0.99636 ), (15.000108899 , 0.00364  ) ],
    'O': [ (15.994914620,  0.99757 ), (16.999131757 , 0.00038  ),
           (17.999159613,  0.00205 )                             ],
    'Si':[ (27.976926535,  0.92223 ), (28.976494665 , 0.04685  ),
           (29.97377001 ,  0.03092 )                             ],
    'P': [ (30.973761998,  1.00    )                             ],
    'S': [ (31.972071174,  0.9499  ), (32.971458910 , 0.0075   ),
           (33.9678670  ,  0.0425  ), (35.967081    , 0.0001   ) ],
    'Na':[ (22.98976928 ,  1.00    )                             ],
    'Cl':[ (34.9688527  ,  0.7576  ), (36.9659026   , 0.2424   ) ],
    'Br':[ (78.918338   ,  0.5069  ), (80.916290    , 0.4931   ) ],
    'Se':[ (73.9224759  ,  0.0089  ), (75.9192137   , 0.0937   ),
           (76.9199142  ,  0.0763  ), (77.917309    , 0.2377   ),
           (79.916522   ,  0.4961  ), (81.916700    , 0.0873   ) ],
    'F': [ (18.99840322 ,  1.0000  ) ],
}
isotopic_distributions_alternative = {
    # proposed alternative - code would need some makeover
    'H': [
        {
            'mass': 1.0078250322,
            'abundance': 0.999885
        },
        {
            'mass': 2.0141017781,
            'abundance': 0.000115
        }
    ]
}

PROTON   = 1.00727646677
ELECTRON = isotopic_distributions['H'][0][0] - PROTON

aa_compositions = {
    # NOTE: This variable will be update if fixed_labels is set.
    #       Ultimately each Isotopoloque class instant will store
    #       the modified dictionary under self.aa_compositions
    'A': 'C3H5NO',
    'C': 'C3H5NOS',
    'D': 'C4H5NO3',
    'E': 'C5H7NO3',
    'F': 'C9H9NO',
    'G': 'C2H3NO',
    'H': 'C6H7N3O',
    'I': 'C6H11NO',
    'K': 'C6H12N2O',
    'L': 'C6H11NO',
    'M': 'C5H9NOS',
    'N': 'C4H6N2O2',
    'P': 'C5H7NO',
    'Q': 'C5H8N2O2',
    'R': 'C6H12N4O',
    'S': 'C3H5NO2',
    'T': 'C4H7NO2',
    'V': 'C5H9NO',
    'W': 'C11H10N2O',
    'Y': 'C9H9NO2',
}
