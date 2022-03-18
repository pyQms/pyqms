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
        * Sarin, L.P.
        * Hippler, M.
        * Leidel, S.A.
        * Fufezan, C.

"""

from __future__ import absolute_import
import sys
import re
import copy
import bisect
import pyqms
import operator
import time
import numpy as np
from chemical_composition import ChemicalComposition


class IsotopologueLibrary(dict):
    """
    The Isotopologue library is the core of pyQms.

    Keyword Arguments:
        molecules (list of str): Molecules used to build the library,
            for more details see below.
        charges (list of int): Charge list used to build the library
        metabolic_labels (dict): see below
        fixed_labels (dict): see below
        params (dict): Match parameters, see `pyqms.params`
        trivial_names (dict): Dictionary that is used to build up lookups.
            Key is a molecule and value a trivial name.
        evidences (dict): Dictionary that is used to build up additional
            lookups. Key is a formula pointing to a subdict. Subdict has
            molecules as keys and values are 'trivial_names' as a list and
            'evidences' holding evidence/identification information
        verbose (bool): Be verbose or not during initialization and matching.

    Keyword argument examples:

    * **molecules** The molecule format can be anything that the
      ChemicalComposition class understands. Currently this can for example
      be::

          [
              '+{0}'.format('H2O'),
              '{peptide}'.format(peptide='PEPTIDE'),
              '{peptide}+{0}'.format('PO3', peptide='PEPTIDE'),
              '{peptide}#{unimod}:{pos}'.format(
                  peptide = peptide,
                  unimod = 'Oxidation',
                  pos = 1
              )
          ]

    * **metabolic_labels** is used to define new element pools
      with enriched isotopes. The dict key defines an enriched element,
      e.g. 15N or 13C and its value is a list of floats [0 - 1.0] defining
      enrichment.The combination of those pools is used to
      calculate isotopologues::

            {
                '15N' : [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
            }

    * **fixed_labels** are based on `unimod`_. Fixed molecules
      do not change the shape of the isotoplogue drastically but introduce
      a simple mass shift, like SILAC, 18O or others.

      The format is for example::

            {
              'R' : ['C(-6) 13C(6) N(-4) 15N(4)','']
            }

    .. _unimod:
        http://www.unimod.org/modifications_list.php?


    Returns:
        dict: Isotopologue library as dict where the top key is always the
        chemical forumla unimod style.

    simplified::

        {
            'C(34)H(53)N(7)O(15)': {
                'cc': {
                    'C': 34, 'H': 53, 'N': 7, 'O': 15
                },
                'env': {
                    (('N', '0.000'),): {
                        # charge :
                        1: {
                            # all transformed mz values
                            'atmzs': {
                                800443,
                                800444,
                                # ... skipped
                                803459,
                                803460,
                                803461
                            },
                            # theoretical mz values
                            'mz': [
                                800.4472772254203,
                                801.450389542063,
                                802.4536114854914,
                                803.4568170275203,
                                804.4597382487398,
                                805.463171867346,
                                806.4631454917885,
                                807.4676949759603
                            ],
                            # transformed mz values within error
                            # packages are on on peak level
                            'tmzs': [
                                {
                                    800443,
                                    800444,
                                    # ... skipped
                                    800451
                                },
                                {
                                    801446,
                                    801447,
                                    # ... skipped
                                    801454
                                },
                                {
                                    802450,
                                    802451,
                                    # ... skipped
                                    802458
                                },
                                {
                                803453,
                                    803454,
                                    803455,
                                    # ... skipped
                                    803461
                                },
                                None,
                                None,
                                None,
                                None
                            ]
                        },
                        # charge independent information
                         'abun': [
                            64799,
                            26251,
                            7164,
                            1456,
                            175,
                            20,
                            1,
                            0
                        ],
                        'c_peak_pos': [
                            0,
                            1,
                            2,
                            3,
                            None,
                            None,
                            None,
                            None
                        ],
                        'isot': [],
                        'mass': [
                            799.3599640346001,
                            800.3629760500413,
                            801.3660976813065,
                            802.3692029128123,
                            803.3720238519379,
                            804.3753571372156,
                            805.3752307742944,
                            806.3796798135622
                        ],
                        'n_c_peaks': 4.0,
                        'relabun': [
                            1.0,
                            0.40511743373159037,
                            0.11054965400744385,
                            0.022466784140529883,
                            0.002693331560888158,
                            0.0003019650321460501,
                            7.716705830708012e-06,
                            3.639837831552297e-08
                        ]
                    }
                }
            }
        }
    """

    def __init__(
        self,
        molecules=None,
        charges=None,
        metabolic_labels=None,
        fixed_labels=None,
        params=None,
        trivial_names=None,
        verbose=True,
        evidences=None,
        unimod_files=None,
        add_default_files=True,
    ):
        assert molecules is not None, "require list of molecules"
        assert charges is not None, "require list of charges"

        self.params = copy.deepcopy(pyqms.params)
        if params is not None:
            self.params.update(params)
        self.zero_labeled_percentile = self.params["PERCENTILE_FORMAT_STRING"].format(0)

        if metabolic_labels is None or metabolic_labels == {}:
            metabolic_labels = {"15N": [0.0]}
        if fixed_labels is None:
            fixed_labels = {}
        self.build_isotoplogues = True
        self.verbose = verbose
        self.metabolic_labels = metabolic_labels
        self.fixed_labels = fixed_labels
        self.charges = charges
        self.lookup = {
            "molecule to formula": {},
            "formula to molecule": {},
            "molecule to annotations": {},
            "molecule fixed label variations": {},
            "formula to trivial name": {},
        }
        self.break_points = {}
        self.skipped_molecules = set()
        if trivial_names is not None:
            self.lookup["molecule to trivial name"] = trivial_names
        if evidences is None:
            evidences = {}
        self.unimod_files = unimod_files
        self.lookup["formula to evidences"] = evidences

        self.regex = {
            "<isotope><element>(<count>)": re.compile(
                r"""
                                                (?P<isotope>[0-9]*)
                                                (?P<element>[A-Z][a-z]*)
                                              \((?P<count>[-]*[0-9]*)\)
                                                """,
                re.VERBOSE,
            ),
            "<isotope><element>": re.compile(
                r"""
                                                (?P<isotope>[0-9]*)
                                                (?P<element>[A-Z][a-z]*)
                                                """,
                re.VERBOSE,
            ),
        }
        self.aa_compositions = {}
        self.isotopic_distributions = {}
        self.computed_level_complex_isotopes = 1
        self.formulas_sorted_by_mz = []
        self.match_sets = {}
        self.match_set_mz_range = [None, None]

        self._cache_kb()  # knowledge_base information
        if self.verbose:
            print("> Metabolic labels        >", self.metabolic_labels)
            print("> Fixed labels            >", self.fixed_labels)
            print("> Charges                 >", self.charges)
            print(
                "> Machine ppm offset      > {MACHINE_OFFSET_IN_PPM}".format(
                    **self.params
                )
            )
        # ----------------------------------------------------------------
        #       METABOLIC LABELS
        # ----------------------------------------------------------------
        self.metabolically_labeled_elements = set()
        # is set in self._extend_isotopic_distributions_with_metabolic_labels()
        self._extend_isotopic_distributions_with_metabolic_labels()
        self.labled_percentiles = []
        self._build_label_percentile_tuples()
        if self.verbose:
            print("> Label percentile tuples >", self.labled_percentiles)
        # ----------------------------------------------------------------
        #       FIXED LABELS
        # ----------------------------------------------------------------
        if len(self.fixed_labels.keys()) != 0:
            # we redefine molecules according to fixed_labeld, i.e SILAC or
            # similar
            self._extend_kb_with_fixed_labels()
            molecules = self._extend_molecules_with_fixed_labels(molecules)

            # will be redefined ...self._highest_element_count
        self._highest_element_count = {}
        self._range_for_elements_with_two_isotopes = [None, None]
        cc_factory = ChemicalComposition(
            aa_compositions=self.aa_compositions,
            isotopic_distributions=self.isotopic_distributions,
            unimod_file_list=unimod_files,
            add_default_files=add_default_files,
        )
        # ----------------------------------------------------------------
        #       BUILDING ISOTOPOLGUES ....
        # ----------------------------------------------------------------
        for molecule in list(set(molecules)):
            if "#" in molecule:
                # molecule is peptide with unimod modification format
                bits = molecule.split("#")
                if len(bits) > 2:
                    raise ValueError(
                        f"{molecule} contains too many '#' {len(bits)} only one allowed"
                    )
                sequence = bits[0]
                modification = bits[1]
            else:
                sequence = molecule
                modification = None
            cc_factory.use(deprecated_format=molecule)
            # mass = cc_factory.mass()
            # if mass / max(self.charges) > self.params['UPPER_MZ_LIMIT']:
            #     continue
            # print(cc_factory.hill_notation_unimod())
            # exit(1)
            chemical_composition = cc_factory.generate_cc_dict()
            # if chemical_composition.mass() / max(self.charges) > self.params['UPPER_MZ_LIMIT']:
            #     self.skipped_molecules.add( molecule )
            #     continue
            # print( chemical_composition.keys() )
            formula = cc_factory.hill_notation_unimod()
            # formula = chemical_composition.hill_notation_unimod()
            if trivial_names is not None:
                molecule_trivial_name = trivial_names.get(molecule, None)
                if molecule_trivial_name is not None:
                    try:
                        self.lookup["formula to trivial name"][formula]
                    except:
                        self.lookup["formula to trivial name"][formula] = []
                    self.lookup["formula to trivial name"][formula].append(
                        molecule_trivial_name
                    )
                else:
                    print("No trivial name for molecule", molecule)
            self.lookup["molecule to formula"][molecule] = formula
            try:
                self.lookup["formula to molecule"][formula]
            except:
                self.lookup["formula to molecule"][formula] = []
            self.lookup["formula to molecule"][formula].append(molecule)
            # this is required since we translate all input molecules
            # (`MAAGALOH+O` or simple `+H2O`) to their respective formula using
            # hill_notation to avoid any mismatches
            if formula not in self.keys():
                self[formula] = {"env": {}, "cc": chemical_composition}
                # self[formula] = {"env": {}, "cc": ChemicalComposition}
            for percentile_tuple in self.labled_percentiles:
                self[formula]["env"][percentile_tuple] = {}

            cc_factory_requires_update = False
            for element, count in chemical_composition.items():
                if element not in self.isotopic_distributions.keys():
                    self.isotopic_distributions[element] = {}
                    # this is a new case where the unimod has new isotopes
                    # that need to be considered ...
                    # see 806
                    pattern = self.regex["<isotope><element>"]
                    match = pattern.match(element)
                    try:
                        enriched_isotope = int(round(float(match.group("isotope"))))
                    except:
                        print("Failed on", element)
                        print("Maybe element is not in pyqms.knowledge_base.py ?")
                        print("Current Distributions available:")
                        print(self.isotopic_distributions.keys())
                        # import csv
                        # ___element_list = []
                        # with open('/home/cf322940/Downloads/NIST_isotope_abundance.csv') as nos:
                        #     for d in csv.DictReader( nos ):
                        #         # print( d.keys())/
                        #         if d['Symbol'] == element:
                        #             ___element_list.append(
                        #                 (
                        #                     d['Relative atomic mass'],
                        #                     d['Relative isotope abundance']
                        #                 )
                        #             )
                        # print('"{0}" : ['.format( element ))
                        # for mass, i in ___element_list:
                        #     print('     ( {0}, {1} ),'.format( mass, i ))
                        # print('],')
                        # print(match)
                        exit(1)
                    # print('> Extending isotopic distribution that are within the modifications (upep)')
                    # enriched_isotope = int(round(float(match.group('isotope'))))
                    template_element = match.group("element")
                    enrichment_key = "{0}{1}".format(enriched_isotope, template_element)
                    target_percentile = (
                        self.params["FIXED_LABEL_ISOTOPE_ENRICHMENT_LEVELS"].get(
                            enrichment_key
                        )
                        or 0.994
                    )

                    print(
                        "> Extending isotopic distribution that are within"
                        " the modifications for "
                        "element {0}, enrichment levels set to {1}".format(
                            enrichment_key, target_percentile
                        )
                    )
                    new_distribution = self._recalc_isotopic_distribution(
                        element=template_element,
                        target_percentile=target_percentile,
                        enriched_isotope=enriched_isotope,
                    )
                    self.isotopic_distributions[element][
                        self.zero_labeled_percentile
                    ] = new_distribution
                    cc_factory_requires_update = True

            if cc_factory_requires_update:
                cc_factory.isotopic_distributions.update(self.isotopic_distributions)
            for element, count in chemical_composition.items():
                try:
                    self._highest_element_count[element]
                except:
                    self._highest_element_count[element] = 0

                if self._highest_element_count[element] < count:
                    self._highest_element_count[element] = count
                for label_percentile in self.isotopic_distributions[element].keys():
                    if len(self.isotopic_distributions[element][label_percentile]) == 2:
                        if (
                            self._range_for_elements_with_two_isotopes[0] is None
                            or count < self._range_for_elements_with_two_isotopes[0]
                        ):
                            self._range_for_elements_with_two_isotopes[0] = count
                        if (
                            self._range_for_elements_with_two_isotopes[1] is None
                            or count > self._range_for_elements_with_two_isotopes[1]
                        ):
                            self._range_for_elements_with_two_isotopes[1] = count
                        break  # we just need to check one :)
                        # No need to iterate over all label_percentiles since they
                        # all share the same number of isotopes
        #
        # here we add 15N count ( or similar ) to natural count in order to
        # facility merge between fixed and metabolic labels
        #
        for element in list(self._highest_element_count.keys()):
            if element.isalpha() is False:
                pattern = self.regex["<isotope><element>"]
                match = pattern.match(element)
                enriched_isotope = int(round(float(match.group("isotope"))))
                template_element = match.group("element")
                if template_element not in self._highest_element_count.keys():
                    self._highest_element_count[template_element] = 0
                self._highest_element_count[
                    template_element
                ] += self._highest_element_count.get(element, 0)
                # print( element, label_percentile, self.isotopic_distributions[ element ])
                for label_percentile in self.isotopic_distributions[
                    template_element
                ].keys():
                    if (
                        len(
                            self.isotopic_distributions[template_element][
                                label_percentile
                            ]
                        )
                        == 2
                    ):
                        # we might need to push the max limit here as well
                        if (
                            self._highest_element_count[template_element]
                            > self._range_for_elements_with_two_isotopes[1]
                        ):
                            self._range_for_elements_with_two_isotopes[
                                1
                            ] = self._highest_element_count[template_element]

            # self._highest_element_count['C'] += 20
            # self._highest_element_count['N'] += 20
        #
        # molecules = molecules - self.skipped_molecules
        # if self.params['verbose']:
        #     print(
        #         'Skipped {0} molecules because their mass (m/z) exceeds the upper m/z limit (m/z limit: {1})'.format(
        #             len(self.skipped_molecules),
        #             self.params['UPPER_MZ_LIMIT']
        #         )
        #     )
        # exit(1)
        self.break_points["Building isotopologues"] = time.time()
        if self.build_isotoplogues:
            self._create_element_trees()
            #
            # Building isotopologues ...
            #
            number_of_formulas = len(self.keys())
            for index, formula in enumerate(self.keys()):
                if self.verbose:
                    print(
                        "> Building isotopologue {0:0>5}/{1:0>5}".format(
                            index + 1, number_of_formulas
                        ),
                        end="\r",
                    )
                # element_stats              = {}
                not_labled_elements = (
                    set(self[formula]["cc"].keys())
                    - self.metabolically_labeled_elements
                )

                pt_dependency = {}
                for label_percentile_tuple in self[formula]["env"].keys():
                    pt_dependency[label_percentile_tuple] = {
                        # this is the not_labeled_but_percentile_dependent_storage
                        "not_labeled_element_combos": [],
                        "zero_isotopic_pos": 0,
                        "element_stats": {},
                        "not_labled_elements": [],
                        "lp_depenent_cc": None,  # is set below
                    }
                    #
                    # here we going to update the element counts
                    # depending on label
                    #
                    lp_dependent_cc = {}
                    for element, count in self[formula]["cc"].items():
                        targe_element_key = element

                        if (
                            element
                            in self.params["FIXED_LABEL_ISOTOPE_ENRICHMENT_LEVELS"]
                        ):
                            # even if not specified in params, self.params
                            # gets updated by _extend_kb_with_fixed_labels
                            isotope_less_element = element[-1]
                            for (
                                percentile_element,
                                label_percentile,
                            ) in label_percentile_tuple:
                                # print(percentile_element, label_percentile, isotope_less_element)
                                if (
                                    percentile_element == isotope_less_element
                                    and label_percentile
                                    == self.params["PERCENTILE_FORMAT_STRING"].format(
                                        self.params[
                                            "FIXED_LABEL_ISOTOPE_ENRICHMENT_LEVELS"
                                        ][element]
                                    )
                                ):
                                    targe_element_key = isotope_less_element
                        # else:
                        #     if element.isalpha() is False:
                        #         continue
                        try:
                            lp_dependent_cc[targe_element_key]
                        except:
                            lp_dependent_cc[targe_element_key] = 0
                        lp_dependent_cc[targe_element_key] += count

                        if targe_element_key not in self.metabolically_labeled_elements:
                            pt_dependency[label_percentile_tuple][
                                "not_labled_elements"
                            ].append(element)

                    pt_dependency[label_percentile_tuple][
                        "lp_depenent_cc"
                    ] = lp_dependent_cc
                    for element, level in lp_dependent_cc.items():
                        if element in self.metabolically_labeled_elements:
                            continue

                        default_element_label = list(
                            self.element_trees[element].keys()
                        )[0]
                        # is this zero_labeled_percentile ?
                        # --------
                        env_range = (
                            self.element_trees[element][default_element_label][level][
                                "maxPos"
                            ]
                            - self.element_trees[element][default_element_label][level][
                                "minPos"
                            ]
                            + 1
                        )
                        pt_dependency[label_percentile_tuple][
                            "not_labeled_element_combos"
                        ].append([element, env_range])
                        pt_dependency[label_percentile_tuple][
                            "zero_isotopic_pos"
                        ] += self.element_trees[element][default_element_label][level][
                            "minPos"
                        ]
                        pt_dependency[label_percentile_tuple]["element_stats"][
                            element
                        ] = (default_element_label, level)

                # -- now the labeled part !! --
                for label_percentile_tuple in self[formula]["env"].keys():

                    current_pt_dep = pt_dependency[label_percentile_tuple]
                    final_local_zero_isotopic_pos = current_pt_dep["zero_isotopic_pos"]
                    labeled_element_combos = []

                    for element, label_percentile in label_percentile_tuple:
                        level = current_pt_dep["lp_depenent_cc"].get(element, 0)
                        if level == 0:
                            continue
                        env_range = (
                            self.element_trees[element][label_percentile][level][
                                "maxPos"
                            ]
                            - self.element_trees[element][label_percentile][level][
                                "minPos"
                            ]
                            + 1
                        )
                        labeled_element_combos.append([element, env_range])
                        current_pt_dep["element_stats"][element] = (
                            label_percentile,
                            level,
                        )
                        final_local_zero_isotopic_pos += self.element_trees[element][
                            label_percentile
                        ][level]["minPos"]

                    tmp = {}
                    for combo in self._create_combinations(
                        current_pt_dep["not_labeled_element_combos"]
                        + labeled_element_combos
                    ):
                        isotope_pos = (
                            sum([pos for element, pos in combo])
                            + final_local_zero_isotopic_pos
                        )
                        # print(combo)
                        try:
                            tmp[isotope_pos]
                        except:
                            tmp[isotope_pos] = {"abun": [], "mass": []}
                        abun = 1
                        mass = 0
                        for element, pos in combo:
                            label_percentile, level = current_pt_dep["element_stats"][
                                element
                            ]
                            min_pos = self.element_trees[element][label_percentile][
                                level
                            ]["minPos"]
                            try:
                                abun *= self.element_trees[element][label_percentile][
                                    level
                                ]["env"][pos + min_pos]["abun"]
                                mass += self.element_trees[element][label_percentile][
                                    level
                                ]["env"][pos + min_pos]["mass"]
                            except:
                                print(
                                    "Can we use the limits ?",
                                    self.element_trees[element][label_percentile][
                                        level
                                    ]["minPos"],
                                    self.element_trees[element][label_percentile][
                                        level
                                    ]["maxPos"],
                                    pos,
                                    element,
                                    combo,
                                    level,
                                    element,
                                    file=sys.stderr,
                                )
                                # print('element',element)
                                # print('self.element_trees[{0}]'.format(element),self.element_trees[element])
                                # print('self.element_trees[{0}][{1}]'.format(element,label_percentile),self.element_trees[element][label_percentile])
                                # print('<>'*40)
                                # print('self.element_trees[{0}][{1}][{2}]'.format(element,label_percentile,level),self.element_trees[element][label_percentile][level])
                                # print( self.element_trees[element][label_percentile][level]['env'][pos] )
                                sys.exit(1)
                        tmp[isotope_pos]["abun"].append(abun)
                        tmp[isotope_pos]["mass"].append(mass)
                    self[formula]["env"][label_percentile_tuple] = {
                        "isot": [],
                        "mass": [],
                        "abun": [],
                        "relabun": [],
                        "c_peak_pos": [],
                        # ^ peaks with higher intensities than
                        # self.params['MIN_REL_PEAK_INTENSITY_FOR_MATCHING']
                        "n_c_peaks": 0,
                        # number of matchable peaks
                    }
                    for charge in self.charges:
                        self[formula]["env"][label_percentile_tuple][charge] = {
                            "mz": [],  # all mz values
                            "tmzs": [],  # transformed mz sets incl. measured precision, pymzml hasPeak style. Each peak one set
                            "atmzs": set(),  # all transformed mz sets together
                        }

                    max_intensity = 0
                    for isotope_pos in sorted(tmp.keys()):
                        total_local_intensity = sum(tmp[isotope_pos]["abun"])
                        if total_local_intensity > max_intensity:
                            max_intensity = total_local_intensity

                    for isotope_pos in sorted(tmp.keys()):
                        total_local_intensity = sum(tmp[isotope_pos]["abun"])
                        if total_local_intensity < sys.float_info.epsilon:
                            continue
                        total_local_mass = 0
                        for n, mass in enumerate(tmp[isotope_pos]["mass"]):
                            total_local_mass += (
                                mass
                                * tmp[isotope_pos]["abun"][n]
                                / float(total_local_intensity)
                            )

                        self[formula]["env"][label_percentile_tuple]["mass"].append(
                            total_local_mass
                        )
                        self[formula]["env"][label_percentile_tuple]["abun"].append(
                            int(
                                round(
                                    total_local_intensity
                                    * self.params["INTENSITY_TRANSFORMATION_FACTOR"]
                                )
                            )
                        )
                        relative_intensity = total_local_intensity / float(
                            max_intensity
                        )
                        self[formula]["env"][label_percentile_tuple]["relabun"].append(
                            relative_intensity
                        )
                        c_peak = False
                        if (
                            relative_intensity
                            >= self.params["MIN_REL_PEAK_INTENSITY_FOR_MATCHING"]
                        ):
                            c_peak = True

                        if c_peak:
                            self[formula]["env"][label_percentile_tuple][
                                "n_c_peaks"
                            ] += 1.0
                            self[formula]["env"][label_percentile_tuple][
                                "c_peak_pos"
                            ].append(isotope_pos)
                            # NOTE: not sure if tuple is needed ...
                            # isotope_pos is not used anywhere after here ?
                        else:
                            self[formula]["env"][label_percentile_tuple][
                                "c_peak_pos"
                            ].append(None)

                        for charge in self.charges:
                            # if charge > 0:
                            #     ionization_spec = pyqms.knowledge_base.PROTON
                            # else:
                            #     ionization_spec = pyqms.knowledge_base.ELECTRON
                            # ^--- negative mode = proton loss not electron addition

                            mz = (
                                total_local_mass + charge * pyqms.knowledge_base.PROTON
                            ) / float(abs(charge))

                            #
                            # MACHINE ERROR
                            #
                            if self.params["MACHINE_OFFSET_IN_PPM"] != 0:
                                mz = (
                                    mz
                                    + mz * 1e-6 * self.params["MACHINE_OFFSET_IN_PPM"]
                                )

                            self[formula]["env"][label_percentile_tuple][charge][
                                "mz"
                            ].append(mz)
                            if c_peak:
                                tmz_set = self._transform_mz_to_set(mz)
                                self[formula]["env"][label_percentile_tuple][charge][
                                    "tmzs"
                                ].append(tmz_set)
                                self[formula]["env"][label_percentile_tuple][charge][
                                    "atmzs"
                                ] |= tmz_set
                            else:
                                self[formula]["env"][label_percentile_tuple][charge][
                                    "tmzs"
                                ].append(None)
                                # now tmzs list has the same length > index n holds

                        #
                        # now add the ranges to the global list
                        #
                    for charge in self.charges:
                        # try:
                        lower_mz = self[formula]["env"][label_percentile_tuple][charge][
                            "mz"
                        ][0]
                        # except:
                        #     print( charge )
                        #     print( formula )
                        #     print( self[ formula ]['env'][label_percentile_tuple] )
                        #     exit(1)
                        lower_mz = self[formula]["env"][label_percentile_tuple][charge][
                            "mz"
                        ][0]
                        upper_mz = self[formula]["env"][label_percentile_tuple][charge][
                            "mz"
                        ][-1]

                        if self.params["LOWER_MZ_LIMIT"] <= lower_mz:
                            if upper_mz <= self.params["UPPER_MZ_LIMIT"]:
                                self.formulas_sorted_by_mz.append(
                                    (
                                        lower_mz,
                                        upper_mz,
                                        charge,
                                        label_percentile_tuple,
                                        formula,
                                    )
                                )
                        # if lower_mz <= self.params['UPPER_MZ_LIMIT']:
                        #     self[ formula ]['env'][label_percentile_tuple]\
                        #         [ charge ] = 'out of mz range'
                    # check if formula is out of range
                    # self.params['UPPER_MZ_LIMIT']
                    # add only peaks that are bigger than
                    # self.params['MOLECULE_MIN_ABUNDANCE']

                # < end > for label_tuple loop

            # < end of molecule look >
        if self.verbose:
            print(
                "\n> Execution time {0:.2f} seconds".format(
                    time.time() - self.break_points["Building isotopologues"]
                )
            )
        # creating formula packages for faster matching
        # requires more ram but is faster
        # alternatively param['MAX_MOLECULES_PER_MATCH_BIN'] can be
        # set to 1
        if self.build_isotoplogues:
            self.formulas_sorted_by_mz.sort()

            number_of_theoretical_formulas = len(self.formulas_sorted_by_mz)
            for package_number, raw_index in enumerate(
                range(
                    0,
                    number_of_theoretical_formulas,
                    self.params["MAX_MOLECULES_PER_MATCH_BIN"],
                )
            ):
                try:
                    next_raw_index = (
                        raw_index + self.params["MAX_MOLECULES_PER_MATCH_BIN"]
                    )
                    self.formulas_sorted_by_mz[next_raw_index]
                except:
                    next_raw_index = number_of_theoretical_formulas
                self.match_sets[package_number] = {
                    "ids": [raw_index, next_raw_index],
                    "tmzs": set(),
                    "mz_range": [None, None],
                }
                for index in range(raw_index, next_raw_index):
                    (
                        lower_mz,
                        upper_mz,
                        charge,
                        label_percentile_tuple,
                        formula,
                    ) = self.formulas_sorted_by_mz[index]
                    if self.params["MAX_MOLECULES_PER_MATCH_BIN"] != 1:
                        self.match_sets[package_number]["tmzs"] |= self[formula]["env"][
                            label_percentile_tuple
                        ][charge]["atmzs"]
                        if (
                            self.match_sets[package_number]["mz_range"][0] is None
                            or lower_mz < self.match_sets[package_number]["mz_range"][0]
                        ):
                            self.match_sets[package_number]["mz_range"][0] = lower_mz
                        if (
                            self.match_sets[package_number]["mz_range"][1] is None
                            or upper_mz > self.match_sets[package_number]["mz_range"][1]
                        ):
                            self.match_sets[package_number]["mz_range"][1] = upper_mz
                    else:
                        self.match_sets[package_number]["tmzs"] = self[formula]["env"][
                            label_percentile_tuple
                        ][charge]["atmzs"]
                        self.match_sets[package_number]["mz_range"][0] = lower_mz
                        self.match_sets[package_number]["mz_range"][1] = upper_mz

                # print(package_number, raw_index, next_raw_index)
                # print(self.formulas_sorted_by_mz[raw_index:next_raw_index])
                # print(self.match_sets[package_number])
            # exit(1)
            # setting global mz_range self.match_set_mz_range
            for package_number in self.match_sets.keys():
                if (
                    self.match_set_mz_range[0] is None
                    or self.match_sets[package_number]["mz_range"][0]
                    < self.match_set_mz_range[0]
                ):
                    self.match_set_mz_range[0] = self.match_sets[package_number][
                        "mz_range"
                    ][0]
                if (
                    self.match_set_mz_range[1] is None
                    or self.match_sets[package_number]["mz_range"][1]
                    > self.match_set_mz_range[1]
                ):
                    self.match_set_mz_range[1] = self.match_sets[package_number][
                        "mz_range"
                    ][1]
            if self.verbose:
                if len(self.formulas_sorted_by_mz) == 0:
                    print("> No molecules have been added into the library ")
                    print("  maybe molecules are out of mz limit range ?")
                    print(
                        "  Current mz limits in params are: [ {LOWER_MZ_LIMIT} .. {UPPER_MZ_LIMIT}]".format(
                            **self.params
                        )
                    )
                    exit(1)
                else:
                    print(
                        "> Created {0} match sets, total mz range [{1:10.5f} .. {2:10.5f}]".format(
                            len(self.match_sets), *self.match_set_mz_range
                        )
                    )
        return

    def _build_label_percentile_tuples(self):
        """
        Builds labled_percentile tuples list, which is stored in
        self.labled_percentiles. Those tuples represent all possible
        combinations of metabolic labels.

        The resulting list has the following format::

            [
                ( ('C','0.000'), ('N','0.000') ),
                ( ('C','0.000'), ('N','0.994') ),
                ( ('C','0.994'), ('N','0.000') ),
                ( ('C','0.994'), ('N','0.994') ),
            ]


        Each tuple in this list contains the combination of the metabolic
        labels.
        """
        pattern = self.regex["<isotope><element>"]
        for combo in self._create_combinations(
            [(k, len(v)) for k, v in self.metabolic_labels.items()]
        ):
            label_tmp_dict = {}
            for entry in combo:
                match = pattern.match(entry[0])
                # isotope = match.group('isotope')
                element = match.group("element")
                index = entry[1]
                label_tmp_dict[element] = self.params[
                    "PERCENTILE_FORMAT_STRING"
                ].format(self.metabolic_labels[entry[0]][index])

            # element_list, label_percentiles = zip(*sorted(label_tmp_dict.items()))
            self.labled_percentiles.append(tuple(sorted(label_tmp_dict.items())))
            # label_percentile_tuple = self._labeled_percentiles_class(**label_tmp_dict)
            # self.label_percentile_tuples.append( label_percentile_tuple )
        return

    def _cache_kb(self):
        """
        Loads default knowledge base information from
        pyqms.knowledge_base.aa_compositions
        pyqms.knowledge_base.isotopic_distributions and stores it under
        self.aa_compositions and self.isotopic_distributions
        """
        default_aa_compositions = pyqms.knowledge_base.aa_compositions
        for aa, composition in default_aa_compositions.items():
            self.aa_compositions[aa] = ChemicalComposition(
                formula="+" + composition, unimod_file_list=self.unimod_files
            )

        for user_aa, composition in self.params.get("AMINO_ACIDS", {}).items():
            self.aa_compositions[user_aa] = ChemicalComposition(
                formula="+" + composition, unimod_file_list=self.unimod_files
            )

        default_isotopic_distributions = pyqms.knowledge_base.isotopic_distributions
        for element, distribution in default_isotopic_distributions.items():
            labeled_percentile = self.zero_labeled_percentile
            self.isotopic_distributions[element] = {labeled_percentile: []}
            for pos, (mass, abundance) in enumerate(distribution):
                self.isotopic_distributions[element][labeled_percentile].append(
                    (mass, abundance, pos)
                )
            self.isotopic_distributions[element][
                labeled_percentile
            ].sort()  # just in case
        return

    def _create_combinations(
        self, input_list, all_combos=None, current_combo=None, pos=0
    ):
        """
        Recursive combination generator.

        Creates a list of all combinations with the respective position of
        an element in a given list. Instead of creating long lists with all
        possible elements, the input just takes the length of the lists and
        returns the position for each list as combinations.

        Examples:

            Input of [("X",2), ("Y",2)] will yield::

                [
                    [("X",0),("Y",0)],
                    [("X",0),("Y",1)],
                    [("X",1),("Y",0)],
                    [("X",1),("Y",1)]
                ]

        Args:
            input_list (list of tuples): First element is the identifier and
                the second element is the length of the list. E.g.
                [("X", number_of_elements_X), ("Y", number_of_elements_Y)]

        Returns:
            all possible combinations of elements of those lists.
        """
        # print( 'cc ',input_list, all_combos)
        if all_combos is None:
            all_combos = []
        if current_combo is None:
            current_combo = []
        if len(input_list) != 0:
            token, number_of_elements = input_list[pos]
            for index in range(number_of_elements):
                extend_list = current_combo + [(token, index)]
                if pos >= len(input_list) - 1:
                    all_combos.append(extend_list)
                else:
                    all_combos = self._create_combinations(
                        input_list,
                        all_combos=all_combos,
                        current_combo=extend_list,
                        pos=pos + 1,
                    )
        return all_combos

    def _create_element_trees(self):
        """
        Creates element trees

        Examples:

           Structure of self.element_trees::

            self.element_trees = {
                '<element>': {
                    '<label incorporation efficiency>': {
                        '<number of atoms>': {
                            'env': {
                                # 0 corresponds to monoisotopic peak
                                # 1 means one additional quasi "neutron"
                                '<isotopic peak position>' : {
                                    'abundance': '<probability of isotopic peak>',
                                    'mass'     : '<mass of isotopic peak>'
                                },
                                # '...' : '...',
                            },
                            'maxPos': '<highest isotopic peak position where peak probability is above params['ELEMENT_MIN_ABUNDANCE']>',
                            'minPos': '<lowest isotopic peak position where peak probability is above params['ELEMENT_MIN_ABUNDANCE']>'
                            },
                        # ...
                        }
                    },
                # ...
                }
        """
        # two elements use binomial distributions and we can cache 'em
        # in fact this nearly always the case ...

        if self._range_for_elements_with_two_isotopes[0] is not None:
            min_number = self._range_for_elements_with_two_isotopes[0]
            max_number = self._range_for_elements_with_two_isotopes[1]
            if self.verbose:
                print(
                    "> Creating binomial cache > [{0} ... {1}]".format(
                        *self._range_for_elements_with_two_isotopes
                    )
                )
            # create factorial_cache
            self._factorial_cache = {0: 1}
            for n in range(1, max_number + 2):
                self._factorial_cache[n] = self._factorial_cache[n - 1] * n

            # create binomial_cache
            self._binomial_cache = {}
            for n in range(min_number, max_number + 1):
                tmpCache = []
                midEnvPos = int(n / 2.0) + 1
                for k in range(midEnvPos):
                    tmpCache.append(
                        self._factorial_cache[n]
                        // (self._factorial_cache[k] * self._factorial_cache[n - k])
                    )
                if n % 2 == 0:  # uneven number of envelope positions
                    self._binomial_cache[n] = (
                        tmpCache + sorted(tmpCache, reverse=True)[1:]
                    )
                else:  # 1 -> even number of envelope positions
                    self._binomial_cache[n] = tmpCache + sorted(tmpCache, reverse=True)
        self.element_trees = {}
        for element, count in self._highest_element_count.items():
            assert (
                element in self.isotopic_distributions.keys()
            ), "Isotopic distribution and abundance of {0} not in knowledge base"
            self.element_trees[element] = {}
            for label_percentile in sorted(
                self.isotopic_distributions[element].keys(), reverse=True
            ):
                self.element_trees[element][label_percentile] = {
                    1: {
                        "maxPos": self.isotopic_distributions[element][
                            label_percentile
                        ][-1][2],
                        "minPos": self.isotopic_distributions[element][
                            label_percentile
                        ][0][2],
                        "env": {},
                    }
                }
                for mass, abundance, envPos in self.isotopic_distributions[element][
                    label_percentile
                ]:
                    self.element_trees[element][label_percentile][1]["env"][envPos] = {
                        "mass": mass,
                        "abun": abundance,
                    }

            for label_percentile in self.isotopic_distributions[element].keys():
                if self.verbose:
                    print(
                        "> Building {0: >5} element tree with a depth of {1:4}, labeling percentile: {2}".format(
                            element, count, label_percentile
                        )
                    )
                number_of_isotopes = len(
                    self.isotopic_distributions[element][label_percentile]
                )
                if number_of_isotopes == 1:
                    for level in range(2, count + 1):
                        previous_level_mass = self.element_trees[element][
                            label_percentile
                        ][level - 1]["env"][0]["mass"]
                        isotope_mass = self.isotopic_distributions[element][
                            label_percentile
                        ][0][0]
                        self.element_trees[element][label_percentile][level] = {
                            "env": {
                                0: {
                                    "abun": 1,
                                    "mass": previous_level_mass + isotope_mass,
                                }
                            },
                            "maxPos": 0,
                            "minPos": 0,
                        }
                elif number_of_isotopes == 2:
                    # power caches, used for mass and abundance
                    # calculation in _increase_element_envelope
                    self._calc_distributions__two_isotopes(
                        element=element, label_percentile=label_percentile, count=count
                    )

                elif number_of_isotopes > 2:
                    self.computed_level_complex_isotopes = 1
                    self._increase_element_envelope(
                        element=element, label_percentile=label_percentile, count=count
                    )
                else:
                    pass
        # import pprint
        # pprint.pprint(self.element_trees['O'])
        # exit(1)
        return

    def _calc_distributions__two_isotopes(
        self, element=None, label_percentile=None, count=None
    ):
        """
        Internal function
        """
        # Hirsch = None
        local_iso_dist = self.isotopic_distributions[element][label_percentile]
        aPowerCache = {1: local_iso_dist[0][1]}
        bPowerCache = {1: local_iso_dist[1][1]}
        aMassPowerCache = {1: local_iso_dist[0][0]}
        bMassPowerCache = {1: local_iso_dist[1][0]}

        beginningZeroK = 0
        min_number_of_elements = self._range_for_elements_with_two_isotopes[0]
        # print( 'calculating two for', element, label_percentile, count )
        for n in range(min_number_of_elements, count + 1):
            # iterate over levels
            self.element_trees[element][label_percentile][n] = {
                "env": {},
                "maxPos": None,
                "minPos": None,
            }
            local_element_tree_pos = self.element_trees[element][label_percentile][n]
            a = local_iso_dist[0][1]
            b = local_iso_dist[1][1]
            aMass = local_iso_dist[0][0]
            bMass = local_iso_dist[1][0]
            somewhereNotZero = False
            endingZero = False
            if n not in self._binomial_cache:
                print(
                    "expected, {0} and got only {1} in element {2}".format(
                        n, max(self._binomial_cache.keys()), element
                    )
                )
                exit(1)

            for k in range(beginningZeroK, len(self._binomial_cache[n])):
                # _binomial_cache[n][k]
                # is 'n choose k'
                # and k = envPos
                local_element_tree_pos["env"][k] = {"mass": None, "abun": 0}
                if not endingZero:
                    try:
                        aPower = aPowerCache[n - k]
                    except:
                        if k != n:
                            try:
                                aPower = aPowerCache[n - k] = a * aPowerCache[n - 1 - k]
                            except:
                                aPower = aPowerCache[n - k] = a ** (n - k)
                        else:
                            aPower = aPowerCache[n - k] = 1
                    try:
                        bPower = bPowerCache[k]
                    except:
                        if k != 0:
                            bPower = bPowerCache[k] = b * bPowerCache[k - 1]
                        else:
                            bPower = bPowerCache[k] = 1
                    try:
                        aMassPower = aMassPowerCache[n - k]
                    except:
                        if k != n:
                            try:
                                aMassPower = aMassPowerCache[n - k] = (
                                    aMass + aMassPowerCache[n - 1 - k]
                                )
                            except:
                                aMassPower = aMassPowerCache[n - k] = aMass * (n - k)
                        else:
                            aMassPower = aMassPowerCache[n - k] = 0
                    try:
                        bMassPower = bMassPowerCache[k]
                    except:
                        if k != 0:
                            bMassPower = bMassPowerCache[k] = (
                                bMass + bMassPowerCache[k - 1]
                            )
                        else:
                            bMassPower = bMassPowerCache[k] = 0
                    local_element_tree_pos["env"][k]["abun"] = (
                        aPower * bPower * self._binomial_cache[n][k]
                    )
                    local_element_tree_pos["env"][k]["mass"] = aMassPower + bMassPower
                    # if element == '(15)N':
                    #     print('Stored mass and abundance for {0}, \
                    #    labeling percentile {1} and lvl {2}'.format( \
                    #    element, label_percentile, n))
                    #     print(local_element_tree_pos['env'])
                    #
                    # Why are all trees build from start to the beginning ?
                    #
                    if (
                        local_element_tree_pos["env"][k]["abun"]
                        <= self.params["ELEMENT_MIN_ABUNDANCE"]
                    ):
                        local_element_tree_pos["env"][k]["abun"] = int(0)
                        if somewhereNotZero:
                            endingZero = True
                        else:
                            beginningZeroK = k + 1
                    else:
                        local_element_tree_pos["maxPos"] = k
                        if not somewhereNotZero:
                            somewhereNotZero = True
                            local_element_tree_pos["minPos"] = k
        return

    def _extend_isotopic_distributions_with_metabolic_labels(self):
        """
        Extends isotopic_distributions with new element pools.

        Metabolic labels introduce a new isotopic distribution.
        This new distribution is calculated using the method
        self._recalc_isotopic_distribution and stored in::

            self.isotopic_distributions['<element>']['<metabolic labeling efficiency>']


        By default labeling percentile 0.000 for N is added. This is the minimum
        required for the isotopologue data structure.


        Format for metabolic_labels is for example::

            {
                '15N' : [0, 0.994],
                '13C' : [0.99]
            }

        """
        for isotope_element, labeled_percentile_list in self.metabolic_labels.items():
            for percentile in labeled_percentile_list:
                if percentile <= sys.float_info.epsilon:
                    continue
                pattern = self.regex["<isotope><element>"]
                match = pattern.match(isotope_element)
                enriched_isotope = int(round(float(match.group("isotope"))))
                template_element = match.group("element")
                new_distribution = self._recalc_isotopic_distribution(
                    element=template_element,
                    target_percentile=percentile,
                    enriched_isotope=enriched_isotope,
                )
                formated_pecentile = self.params["PERCENTILE_FORMAT_STRING"].format(
                    percentile
                )
                self.isotopic_distributions[template_element][
                    formated_pecentile
                ] = new_distribution
                self.metabolically_labeled_elements.add(template_element)

        if len(self.metabolically_labeled_elements) == 0:
            # nothing was added - no metabolic labeling specified
            # then we add our default zero label element
            self.metabolically_labeled_elements.add("N")

    def _extend_kb_with_fixed_labels(self, synthetic_abundance=0.994):
        """
        Updates knowledge base dictionaries with fixed_label properties.

        The knowledge base is extended with newly defined amino acids and
        elements showing synthetic isotopic distributions. By default
        synthetic_abundance is taken from
        params['FIXED_LABEL_ISOTOPE_ENRICHMENT_LEVELS'] and **fall back**
        value is an abundance of 0.994. This is done using the
        function self._recalc_isotopic_distribution

        format for fixed_labels is for example::

            {
                'R' : ['C(-6) 13C(6) N(-4) 15N(4)',''],
                'K' : '...'
            }

        Heavy Arg will become R0 if defined at index 0,
        i.e. self.fixed_labels['R'][0]

        Note:
            self.aa_compositions and self.isotopic_distributions
            are based on pyqms.knowledge_base and cached by `self._cache_kb`
        """
        for labeled_aa, label_definitions in self.fixed_labels.items():
            for pos, uni_mod_string in enumerate(label_definitions):
                new_aa = "{0}{1}".format(labeled_aa, pos)
                try:
                    new_aa_composition = copy.deepcopy(self.aa_compositions[labeled_aa])
                except:
                    print("Error in _extend_kb_with_fixed_labels")
                    print(self.aa_compositions)
                    print(labeled_aa)
                    # print( self.aa_compositions[ labeled_aa ] )
                    # print( new_aa_composition = copy.deepcopy(
                    #     self.aa_compositions[ labeled_aa ] )
                    # )
                    exit(1)
                formated_umod_list = self.regex["<isotope><element>(<count>)"].findall(
                    uni_mod_string
                )
                for isotope, element, count in formated_umod_list:
                    if isotope == "":
                        formated_element = element
                        # labeling is metabolically added,
                        # i.e. follows pulse chase if needed
                    else:
                        isotope_sorted_by_abundance = sorted(
                            self.isotopic_distributions[element][
                                self.zero_labeled_percentile
                            ],
                            key=operator.itemgetter(1),
                            reverse=True,
                        )
                        natural_isotope = isotope_sorted_by_abundance[0]
                        natural_isoptop_trivial_name = str(
                            round(natural_isotope[0])
                        ).split(".")[0]
                        natural_isotope_percentile = natural_isotope[1]
                        if isotope == natural_isoptop_trivial_name:
                            # this is chemically added and should not follow
                            # pulse chase ...
                            percentile = 0.0  # no enrichment .. stay as you are ..
                            isotope = natural_isoptop_trivial_name

                        # print( self.metabolic_labels )
                        # print( str(mass_of_first_isotope).split('.')[0] )
                        # print( new_aa_composition )
                        # print( isotope_sorted_by_abundance )
                        # formated_element = element
                        # if
                        else:
                            percentile = synthetic_abundance
                        #
                        formated_element = "{0}{1}".format(isotope, element)
                        # print(percentile)
                        if (
                            formated_element
                            in self.params[
                                "FIXED_LABEL_ISOTOPE_ENRICHMENT_LEVELS"
                            ].keys()
                        ):
                            percentile = self.params[
                                "FIXED_LABEL_ISOTOPE_ENRICHMENT_LEVELS"
                            ][formated_element]
                        else:
                            if self.verbose:
                                print(
                                    ">Enrichment levels of {0} was not specified in params."
                                    "Falling back to {1} and inserting temporarily into self.params".format(
                                        formated_element, percentile
                                    )
                                )
                            self.params["FIXED_LABEL_ISOTOPE_ENRICHMENT_LEVELS"][
                                formated_element
                            ] = percentile

                    new_aa_composition[formated_element] += int(count)
                    if isotope != "":
                        enriched_isotope = isotope
                        template_element = element
                        # percentile       = synthetic_abundance
                        new_distribution = self._recalc_isotopic_distribution(
                            element=template_element,
                            target_percentile=percentile,
                            enriched_isotope=enriched_isotope,
                        )
                        formated_pecentile = self.params[
                            "PERCENTILE_FORMAT_STRING"
                        ].format(percentile)
                        if formated_element not in self.isotopic_distributions.keys():
                            self.isotopic_distributions[formated_element] = {}
                        self.isotopic_distributions[formated_element][
                            formated_pecentile
                        ] = new_distribution
                to_be_deleted = []
                for k, v in new_aa_composition.items():
                    if v == 0:
                        to_be_deleted.append(k)
                for k in to_be_deleted:
                    del new_aa_composition[k]
                # print(new_aa_composition)
                self.aa_compositions[new_aa] = new_aa_composition
        return

    def _extend_molecules_with_fixed_labels(self, molecules):
        # self.params['SILAC_AAS_LOCKED_IN_EXPERIMENT']  = ['K','R']
        """
        Reformat molecules depending on the fixed_label.

        All possible combinations of the fixed labels are produced and the
        Molecule names are reformated.

        Examples:
            fixed_labels input::

                {
                    'R': ['C(-6) 13C(6)', ''],
                    # Arg in 13C (pos 0) and without label (pos1)
                    'C': ['H(4) C(2) O S']
                    # DeStreak or Cysteine mercaptoethanol
                }

            Given the peptide EILCEWRRAR, this will produce the following 8
            peptides:

                #. EILC\ **0**\ EWR\ **0**\ R\ **0**\ AR\ **0**
                #. EILC\ **0**\ EWR\ **1**\ R\ **0**\ AR\ **0**
                #. EILC\ **0**\ EWR\ **0**\ R\ **1**\ AR\ **0**
                #. EILC\ **0**\ EWR\ **1**\ R\ **1**\ AR\ **0**
                #. EILC\ **0**\ EWR\ **0**\ R\ **0**\ AR\ **1**
                #. EILC\ **0**\ EWR\ **1**\ R\ **0**\ AR\ **1**
                #. EILC\ **0**\ EWR\ **0**\ R\ **1**\ AR\ **1**
                #. EILC\ **0**\ EWR\ **1**\ R\ **1**\ AR\ **1**

        Note:
            The param entry SILAC_AAS_LOCKED_IN_EXPERIMENT allows amino acids
            to be locked with a label in certain combinations, so that they
            will always appear together.
        """
        # [('9R', 0), ('7R', 0), ('6R', 0), ('3C', 0)] EILC0EWR0R0AR0 EILCEWRRAR
        # [('9R', 0), ('7R', 0), ('6R', 1), ('3C', 0)] EILC0EWR1R0AR0 EILCEWRRAR
        # [('9R', 0), ('7R', 1), ('6R', 0), ('3C', 0)] EILC0EWR0R1AR0 EILCEWRRAR
        # [('9R', 0), ('7R', 1), ('6R', 1), ('3C', 0)] EILC0EWR1R1AR0 EILCEWRRAR
        # [('9R', 1), ('7R', 0), ('6R', 0), ('3C', 0)] EILC0EWR0R0AR1 EILCEWRRAR
        # [('9R', 1), ('7R', 0), ('6R', 1), ('3C', 0)] EILC0EWR1R0AR1 EILCEWRRAR
        # [('9R', 1), ('7R', 1), ('6R', 0), ('3C', 0)] EILC0EWR0R1AR1 EILCEWRRAR
        # [('9R', 1), ('7R', 1), ('6R', 1), ('3C', 0)] EILC0EWR1R1AR1 EILCEWRRAR
        # Found R at 6 R6:['C(-6) 13C(6)', '']
        # Found R at 7 R7:['C(-6) 13C(6)', '']
        # Found R at 9 R9:['C(-6) 13C(6)', '']
        # Found C at 3 C3:['C(1)O(2)']
        #
        # aas_locked_in_exp should go into params !!
        # Maybe have a nicer name ?
        #
        extended_set_of_molecules = set()
        patterns = {}
        for aa in self.fixed_labels.keys():
            patterns[aa] = re.compile(aa)
        if self.params["SILAC_AAS_LOCKED_IN_EXPERIMENT"] is not None:
            assert (
                len(self.params["SILAC_AAS_LOCKED_IN_EXPERIMENT"]) is not 1
            ), 'only one aa in self.params["SILAC_AAS_LOCKED_IN_EXPERIMENT"]  ?'
            all_locked_aa_fixed_mods_have_same_length = True
            for n in range(1, len(self.params["SILAC_AAS_LOCKED_IN_EXPERIMENT"])):
                if len(
                    self.fixed_labels[self.params["SILAC_AAS_LOCKED_IN_EXPERIMENT"][0]]
                ) != len(
                    self.fixed_labels[self.params["SILAC_AAS_LOCKED_IN_EXPERIMENT"][n]]
                ):
                    all_locked_aa_fixed_mods_have_same_length = False
            assert (
                all_locked_aa_fixed_mods_have_same_length is True
            ), "Length of fixed_labels for aas that ought to be locked into \
                same experiment dont show same number of labels ..."

        for full_molecule in molecules:
            label_combinations = []
            haz_mod = False
            for peptide_en_pos, aa in enumerate(full_molecule):
                if aa not in self.aa_compositions.keys():
                    haz_mod = True
                    break
            if haz_mod:
                molecule = full_molecule[:peptide_en_pos]
                addon = full_molecule[peptide_en_pos:]
            else:
                molecule = full_molecule
                addon = ""
            # print('>>>',full_molecule, molecule, peptide_en_pos,addon)
            for aa, pattern in patterns.items():
                for match in pattern.finditer(molecule):
                    label_combinations.append(
                        ("{0}{1}".format(match.start(), aa), len(self.fixed_labels[aa]))
                    )
                    # print('>>>>>>>>>>>>>>>>Found',aa,'at',match.start(), '{0}{1}:{2}'.format(aa,match.start(),self.fixed_labels[aa]))
                    # exit(1)
            if len(label_combinations) == 0:
                extended_set_of_molecules.add(molecule + addon)
            else:
                label_combinations = self._create_combinations(
                    sorted(label_combinations, reverse=True)
                )
                for combination in label_combinations:
                    # print('Combo:',combination)
                    modified_molecule = molecule
                    states = set()
                    sorted_exchange_list = []
                    for element in combination:
                        index_in_molecule = int(element[0][:-1])
                        aa = element[0][-1]
                        fixed_label_index = element[1]
                        sorted_exchange_list.append(
                            (index_in_molecule, aa, fixed_label_index)
                        )
                    sorted_exchange_list.sort(reverse=True)
                    for (
                        index_in_molecule,
                        aa,
                        fixed_label_index,
                    ) in sorted_exchange_list:
                        modified_molecule = (
                            "{seq1}{aa}{fixed_label_index}{seq2}".format(
                                seq1=modified_molecule[:index_in_molecule],
                                seq2=modified_molecule[index_in_molecule + 1 :],
                                aa=aa,
                                fixed_label_index=fixed_label_index,
                            )
                        )
                        if (
                            self.params["SILAC_AAS_LOCKED_IN_EXPERIMENT"] is not None
                            and aa in self.params["SILAC_AAS_LOCKED_IN_EXPERIMENT"]
                        ):
                            states.add(fixed_label_index)
                        #     print(modified_molecule, 'aa', aa, 'in SILAC_LOCK')
                        # else:
                        #     print(modified_molecule, 'not added' )
                    extend_variation_lookup = False
                    modified_molecule_incl_addon = modified_molecule + addon
                    if self.params["SILAC_AAS_LOCKED_IN_EXPERIMENT"] is not None:
                        if len(states) == 1:
                            all_states_in_molecule = set()
                            state_pattern = re.compile(
                                r"""(?P<AA>[A-Z]{1})(?P<state>[0-9]+)""", re.VERBOSE
                            )
                            for mol_state_aa, mol_state_pos in state_pattern.findall(
                                modified_molecule
                            ):
                                if (
                                    mol_state_aa
                                    in self.params["SILAC_AAS_LOCKED_IN_EXPERIMENT"]
                                ):
                                    all_states_in_molecule.add(int(mol_state_pos))
                            if all_states_in_molecule == states:
                                extended_set_of_molecules.add(
                                    modified_molecule_incl_addon
                                )
                                extend_variation_lookup = True
                    else:
                        extended_set_of_molecules.add(modified_molecule_incl_addon)
                        extend_variation_lookup = True
                    if extend_variation_lookup:
                        try:
                            self.lookup["molecule fixed label variations"][
                                full_molecule
                            ]
                        except:
                            self.lookup["molecule fixed label variations"][
                                full_molecule
                            ] = set()
                        self.lookup["molecule fixed label variations"][
                            full_molecule
                        ].add(modified_molecule_incl_addon)
                    # print('>>>',
                    #     molecule,
                    #     modified_molecule,
                    #     index_in_molecule,
                    #     modified_molecule[:index_in_molecule],
                    #     aa
                    # )
        return extended_set_of_molecules

    def _increase_element_envelope(self, element=None, label_percentile=None, count=0):
        """
        Calculates envelopes of an element that has more than two isotopes.

        Args:
            element (str): element to consider (e.g. Oxygen)
            count (int): maximum tree depth i.e. maximum element count

        Stores element tree in self.element_trees
        """
        # print("Debug _increase_element_envelope with element/count", element, count)
        if count != 0:
            count -= 1
            # for  label_percentile  in self.isotopic_distributions[ element ].keys():
            n = self.computed_level_complex_isotopes + 1
            # import pprint
            # pprint.pprint(self.element_trees[element])
            # if  label_percentile  not in self.element_trees[element]:
            #     self.element_trees[element][ label_percentile ] = {}
            self.element_trees[element][label_percentile][n] = {
                "env": {},
                "minPos": None,
                "maxPos": None,
            }
            # print('>>>',n)
            # print('FIX THE BUG ....')
            # print( self.element_trees[element][ label_percentile ] )
            # print( "self.isotopic_distributions[element]", self.isotopic_distributions[element])
            # print( "element", element )
            # print( " self.element_trees[element] ", self.element_trees[element] )
            for envPos in sorted(
                self.element_trees[element][label_percentile][
                    self.computed_level_complex_isotopes
                ]["env"].keys()
            ):
                for (
                    isotopeMass,
                    isotopeAbundance,
                    isotopePos,
                ) in self.isotopic_distributions[element][label_percentile]:
                    if (
                        envPos + isotopePos
                        not in self.element_trees[element][label_percentile][n][
                            "env"
                        ].keys()
                    ):
                        self.element_trees[element][label_percentile][n]["env"][
                            envPos + isotopePos
                        ] = {"mass": [], "abun": 0}

                    # abundance
                    previousAbundance = self.element_trees[element][label_percentile][
                        self.computed_level_complex_isotopes
                    ]["env"][envPos]["abun"]
                    isotopeAbundance = self.isotopic_distributions[element][
                        label_percentile
                    ][isotopePos][1]
                    self.element_trees[element][label_percentile][n]["env"][
                        envPos + isotopePos
                    ]["abun"] += (previousAbundance * isotopeAbundance)

                    # mass
                    previousMass = self.element_trees[element][label_percentile][
                        self.computed_level_complex_isotopes
                    ]["env"][envPos]["mass"]
                    isotopeMass = self.isotopic_distributions[element][
                        label_percentile
                    ][isotopePos][0]
                    self.element_trees[element][label_percentile][n]["env"][
                        envPos + isotopePos
                    ]["mass"].append(previousMass + isotopeMass)

            minPos = False
            for envPos in self.element_trees[element][label_percentile][n][
                "env"
            ].keys():  # calc mean mass
                pathsMassSum = sum(
                    self.element_trees[element][label_percentile][n]["env"][envPos][
                        "mass"
                    ]
                )
                pathsNumber = len(
                    self.element_trees[element][label_percentile][n]["env"][envPos][
                        "mass"
                    ]
                )
                self.element_trees[element][label_percentile][n]["env"][envPos][
                    "mass"
                ] = float(pathsMassSum) / float(pathsNumber)
                if not minPos:
                    if (
                        self.element_trees[element][label_percentile][n]["env"][envPos][
                            "abun"
                        ]
                        > self.params["ELEMENT_MIN_ABUNDANCE"]
                    ):
                        # this will find lowest and highest env pos where abundance is > 0
                        self.element_trees[element][label_percentile][n][
                            "minPos"
                        ] = envPos  #  the range between these positions is ignored
                        minPos = True
                if (
                    self.element_trees[element][label_percentile][n]["env"][envPos][
                        "abun"
                    ]
                    > self.params["ELEMENT_MIN_ABUNDANCE"]
                ):
                    self.element_trees[element][label_percentile][n]["maxPos"] = envPos
            self.computed_level_complex_isotopes += 1
            self._increase_element_envelope(
                element=element, count=count, label_percentile=label_percentile
            )
            # efficiencyList = efficiencyList)
        return

    def match_all(
        self, mz_i_list=None, file_name=None, spec_id=None, spec_rt=None, results=None
    ):
        """
        Matches all isotopologues in the library agains a given *mz_i_list*

        Args:
            mz_i_list (list of tuples): Spectrum information that should be
                matched against. Tuples of m/z and intensity
            file_name (str): Information used for storage purpose. Useful if
                multiple files are parsed with one pyqms.result instance.
            spec_id (int): Information used for storage purpose.
            spec_rt (float): Information used for storage purpose.
            results (`pyqms.Results`): (optional)


        If a results object is passed to match_all, then this object will be
        updated and returned. This is for e.g. to accumulate results for a whole
        LC-MS/MS run.

        For various examples using match_all please refer to the example scripts.

        Returns:

            results class object (obj): Object holding all quantitative information

        """
        if results is None:
            results = pyqms.Results(
                aa_compositions=self.aa_compositions,
                charges=self.charges,
                fixed_labels=self.fixed_labels,
                isotopic_distributions=self.isotopic_distributions,
                lookup=self.lookup,
                metabolic_labels=self.metabolic_labels,
                params=self.params,
            )
        lower_value = (self.match_set_mz_range[0], 0)
        upper_value = (self.match_set_mz_range[1], 0)
        borders = (lower_value, upper_value)
        sliced_spec = self._slice_list(mz_i_list, borders)
        for package_number in self.match_sets.keys():
            mz_range = self.match_sets[package_number]["mz_range"]
            spec_tmz_set, spec_tmz_lookup = self._transform_spectrum(
                sliced_spec, mz_range=mz_range
            )
            # print(self.match_sets[package_number]['tmzs'])
            # print(spec_tmz_set)
            if (
                len(spec_tmz_set & self.match_sets[package_number]["tmzs"])
                >= self.params["MINIMUM_NUMBER_OF_MATCHED_ISOTOPOLOGUES"]
            ):
                for index in range(
                    self.match_sets[package_number]["ids"][0],
                    self.match_sets[package_number]["ids"][1],
                ):
                    # >>>
                    match_results = self.match_isotopologue(
                        index=index,
                        spec_tmz_set=spec_tmz_set,
                        spec_tmz_lookup=spec_tmz_lookup,
                        mz_score_percentile=results.params["MZ_SCORE_PERCENTILE"],
                    )
                    if match_results is None:
                        continue
                    score, scaling_factor, matched_peaks = match_results
                    if score < self.params["M_SCORE_THRESHOLD"]:
                        continue
                    number_of_matched_peaks = 0
                    for mmz, mi, ri, cmz, ci in matched_peaks:
                        if mmz is not None:
                            number_of_matched_peaks += 1
                    if (
                        number_of_matched_peaks
                        < self.params["MINIMUM_NUMBER_OF_MATCHED_ISOTOPOLOGUES"]
                    ):
                        continue

                    (
                        lower_mz,
                        upper_mz,
                        charge,
                        label_percentile_tuple,
                        formula,
                    ) = self.formulas_sorted_by_mz[index]
                    key = (file_name, formula, charge, label_percentile_tuple)
                    value = (spec_id, spec_rt, score, scaling_factor, matched_peaks)
                    # print('added', score, matched_peaks )
                    results.add(key, value)
        return results

    def match_isotopologue(
        self,
        index=None,
        formula=None,
        charge=None,
        label_percentile=None,
        spec_tmz_set=None,
        spec_tmz_lookup=None,
        mz_i_list=None,
        mz_score_percentile=None,
    ):
        """
        Matches a single isotopologue onto a *mz_i_list* or *spec_tmz_set*

        Args:
            index (int): Using this index one can retrieve all information
                about the molecule, i.e. lower_mz, upper_mz, charge,
                label_percentile, formula from *self.formulas_sorted_by_mz*.
                Alternatively, one can use the more
                verbose option: formula, charge and label_percentile
            formula (str): pyqms formula type
            charge (int): molecule charge
            label_percentile: pyQms label percentile
            mz_i_list (list of tuples): List of m/z and intensity tuples, will
                be transformed to a spec_tmz_set given the defined precession.
                Alternatively, spec_tmz_set can be used as input.
            spec_tmz_set (set of ints): tmz value set used for matching.
                Requires spec_tmz_lookup to get the actual mz which is required
                for scoring.
            mz_score_percentile (float): Weighting of mz used for scoring.
                (1 - mz_score_percentile) is then intensity weighting.
                Values 0 - 1.0.

        Note:
            Depending on the machine (some measure intensity better than others)
            adjusting mz_score_percentile value will give
            more accurate results. Best adjusted in pyqms.params (which
            can be passed during isotoplogue lib initialization)

        Returns:
            Match results (tuple of score, scaling factor and matched peaks).

                * **score** reflects the fit of the theoretical isotopologue to the
                  measured (both mz and intensities are compared)
                * **scaling factor** reflects the actual amount of the molecule in the
                  respective spectrum. It is defined as the sum of the total measured
                  intensities divided by the sum of the total calculated intensities
                * **matched_peaks** is list of tuples that contain
                  measured_mz, measured_i, rel_i, calculated_mz, calculated_i

        Multiple m/z values can occur in the range of the measured precision
        of every peak of the isotopologue, thus all combinations are
        considered and scored. Only the best scored match is returned for
        each isotopologue.

        """
        if index is None:
            assert formula is not None, "require formula information for match"
            assert charge is not None, "require charge information for match"
            assert label_percentile is not None, "require tuple list"
            lower_mz = self[formula]["env"][label_percentile][charge]["mz"][0]
            upper_mz = self[formula]["env"][label_percentile][charge]["mz"][-1]
            # really easier via the index :)
            # use function blabla (to be written) to scan
            # self.formulas_sorted_by_mz for your target(s)
        else:
            (
                lower_mz,
                upper_mz,
                charge,
                label_percentile,
                formula,
            ) = self.formulas_sorted_by_mz[index]
        if mz_i_list is not None:
            assert (
                spec_tmz_set is None
            ), "both a transformed spec and mz_i_list were given ... which one to use ?"
            borders = ((lower_mz, 0), (upper_mz, 0))
            sliced_spec = self._slice_list(mz_i_list, borders)
            if spec_tmz_set is None:
                spec_tmz_set, spec_tmz_lookup = self._transform_spectrum(
                    sliced_spec, mz_range=(lower_mz, upper_mz)
                )
        full_isotopologue_set = self[formula]["env"][label_percentile][charge]["atmzs"]
        overlap = full_isotopologue_set & spec_tmz_set
        n_c_peaks = self[formula]["env"][label_percentile]["n_c_peaks"]
        results = None
        match_it = True
        if len(overlap) < self.params["MINIMUM_NUMBER_OF_MATCHED_ISOTOPOLOGUES"]:
            match_it = False
        if len(overlap) / n_c_peaks < self.params["REQUIRED_PERCENTILE_PEAK_OVERLAP"]:
            match_it = False
        if match_it:
            # print('found {0} % matches'.format( len(overlap)/n_c_peaks))
            # print('Overlap:',overlap)
            # NOTE: multiple peaks might match in a given error the theoretical
            # peak, i.e. spec_tmz_lookup contains multiple value.
            #
            # Therefore all possible combinations of those hits have to
            # be evaluate in order to maximize the fit
            matched_mmz_on_isotope_pos = {}
            # sorted_overlap = sorted(  overlap )
            matched_peaks = {}

            for n, isotope_pos in enumerate(
                self[formula]["env"][label_percentile]["c_peak_pos"]
            ):
                if isotope_pos is None:
                    continue
                matched_peaks[n] = [
                    None,
                    None,
                    self[formula]["env"][label_percentile]["relabun"][n],
                    self[formula]["env"][label_percentile][charge]["mz"][n],
                    self[formula]["env"][label_percentile]["abun"][n],
                ]
                tmz = self[formula]["env"][label_percentile][charge]["tmzs"][n]
                matched_mmz_on_isotope_pos[n] = list(overlap & tmz)

            match_combinations = []
            for n, match_list in sorted(matched_mmz_on_isotope_pos.items()):
                if len(match_list) > 0:
                    match_combinations.append((n, len(match_list)))
            # print('> Match combos:', match_combinations )
            scores = []
            for match_combo in self._create_combinations(match_combinations):
                # print( match_combo )
                for isotope_pos, match_index in match_combo:
                    matched_tmz = matched_mmz_on_isotope_pos[isotope_pos][match_index]
                    measured_mz = spec_tmz_lookup[matched_tmz][0][0]
                    measured_i = spec_tmz_lookup[matched_tmz][0][1]
                    matched_peaks[isotope_pos][0] = measured_mz
                    matched_peaks[isotope_pos][1] = measured_i
                # print( matched_peaks.values() , match_combo )
                score, scaling_factor = self.score_matches(
                    matched_peaks.values(), mz_score_percentile
                )
                scores.append(
                    (
                        score,
                        scaling_factor,
                        tuple(
                            (mmz, mi, ri, cmz, ci)
                            for mmz, mi, ri, cmz, ci in matched_peaks.values()
                        ),
                    )
                )
            scores.sort(reverse=True)
            try:
                results = scores[0]
            except:
                results = None
            # print( results )
            # import pprint
            # pprint.pprint([ (so, sc, mp) for so, sc, mp in  scores ])
        return results

    def print_overview(self, formula, charge=None):
        """
        Prints an overview of a given molecule or formula to the std.out

        Args:
            formula (str): Either formula or molecule
            charge (int): Charge of the molecule

        Examples:
            For PEPTIDE and charge 1::

                Chemical formula C(34)H(53)N(7)O(15)
                (('N', '0.000'),)
                Isotope                                               Abundance
                pos      Mass          m/z [MH]+1               transformed    rel.
                 0   799.3599640346      800.4472772254         64799  1.00000000000   0
                 1   800.3629760500      801.4503895421         26251  0.40511743373   1
                 2   801.3660976813      802.4536114855          7164  0.11054965401   2
                 3   802.3692029128      803.4568170275          1456  0.02246678414   3
                 4   803.3720238519      804.4597382487           175  0.00269333156   None
                 5   804.3753571372      805.4631718673            20  0.00030196503   None
                 6   805.3752307743      806.4631454918             1  0.00000771671   None
                 7   806.3796798136      807.4676949760             0  0.00000003640   None


        """
        import pprint

        if formula not in self.keys():
            # maybe we got original molecule as input
            formula = self.lookup["molecule to formula"][formula]
        # assert for
        if charge is None:
            assert isinstance(charge, int), "charge has to be an integer"
            charge = self.charges[0]
        print("\n\n{0: ^30}\n\n".format("Overview"))
        print("> Chemical formula", formula)
        tl = self.lookup["formula to trivial name"].get(formula, None)
        if tl is not None:
            print("> Trivial name{0} {1}".format("" if len(tl) == 1 else "s", tl))
        for label_percentile in sorted(self[formula]["env"]):
            print("> Label percentile", label_percentile)
            print(
                """> Isotope pattern                                                Abundance\n pos      Mass\t\t\tm/z [MH]{0:+1}               transformed    rel.""".format(
                    charge
                )
            )
            # max_intensity = max(self[ formula ]['env'][ label_percentile ]['abun'])
            number_of_isotopic_peaks = len(
                self[formula]["env"][label_percentile]["abun"]
            )
            for n in range(number_of_isotopic_peaks):
                mass = self[formula]["env"][label_percentile]["mass"][n]
                abun = self[formula]["env"][label_percentile]["abun"][n]
                mzc1 = self[formula]["env"][label_percentile][charge]["mz"][n]
                relabun = self[formula]["env"][label_percentile]["relabun"][n]
                c_pos = self[formula]["env"][label_percentile]["c_peak_pos"][n]
                print(
                    "{0: >3}\t{1:16.10f}\t{2:16.10f}\t{3:10.0f}\t{4:12.11f}\t{5}".format(
                        n, mass, mzc1, abun, relabun, c_pos
                    )
                )

        # return
        # pprint.pprint( self[ formula ]['env'][ label_percentile ]['c_peak_pos'] )
        # exit(1)

    def _recalc_isotopic_distribution(
        self, element=None, target_percentile=None, enriched_isotope=None
    ):
        """
        Calculate new isotopic distribution for an enriched element.

        Natural abundance is scaled by (1 - target_percentile), where as
        the target_percentile is added onto the enriched_isotope of choice.
        Ensures that the other abundances are adjusted according to the
        remaining abundance (e.g. 99.4 percent Nitrogen enrichment, 0.6 percent
        remaining natural abundance).

        Args:
            element (str): Element that is used as template.
            target_percentile (float): Enrichment level [0 - 1.0]
            enriched_isotope (int): enriched_isotope is int of mass,
                e.g. 13 for '13C'

        Returns:
            new_distribution (list): List containing tuples of mass, abundance
                and peak position.

            * **mass** original mass of template element
            * **abundance** recalculated abundance
            * **peak position** peak position


        """
        # print('recalc', element, target_percentile, enriched_isotope )
        new_distribution = []
        # target_percentile -=

        total_other_isotope_abundance = 0

        for mass, abundance, pos in self.isotopic_distributions[element][
            self.zero_labeled_percentile
        ]:
            if int(round(mass)) != int(enriched_isotope):
                total_other_isotope_abundance += abundance
                continue
            natural_abundance = abundance

        diff_in_abundance = natural_abundance - target_percentile
        for mass, abundance, pos in self.isotopic_distributions[element][
            self.zero_labeled_percentile
        ]:

            share_in_difference = (
                abundance * diff_in_abundance / total_other_isotope_abundance
            )

            if int(round(mass)) == int(enriched_isotope):
                abundance = target_percentile
            else:
                abundance += share_in_difference
            new_distribution.append((mass, abundance, pos))
        return new_distribution

    def score_matches(self, matched_peaks, mz_score_percentile):
        """
        Score matched peaks.

        Args:
            matched_peaks (list of tuples): List of tuples containing

                * measured_mz (mmz)
                * measured_intensity (mi)
                * relative_intensity_of_calculated_isotopologue_peak (ri)
                * calculated_mz (cmz)
                * calculated_i (ci)

            mz_score_percentile (float): weighting of mz score

        Parameters that influence the scoring are 'MIN_REL_PEAK_INTENSITY_FOR_MATCHING'

        .. _Gower:
            http://venus.unive.it/romanaz/modstat_ba/gowdis.pdf



        **Example plots**

        The figure below highlights the scoring principle. Erros for m/z and
        intensity values are determined and combined into the final mScore. For
        each peak of the isotopologue both errors are determined and influence
        the final score.


        .. image:: images/scoring_principle.png


        Calculated intensities are scaled to match the measured value and the
        deviation is calculated. The lower the intensity, the less accurate
        teh actual peaks are represented. To compensate for this, the intensity
        score decreases faster for large relative intensities compared to small
        relative intensities. This is highlighted in the following figure.
        Legend, x-axis represents the relative intensity error (measured -
        theoretical intensity) and the y-axis the intensity score. Different
        colors represent various relative peak intensities.


        .. image:: images/Intensity_score.png



        **Scoring**

        .. note::

            The proper display of the formulas of the next section
            requires access to the Internet when browsing the HTML
            documentation. The formulas are correctly embedded into the pdf of
            the documentation.


        The pyQms matching score (mScore) is based on the work of  `Gower`_ (1971)
        *A General Coefficient of Similarity and Some of Its Properties*,
        Biometrics (27), 857-871.
        The matching and scoring is performed on the m/z values and the
        intensity values independently yielding two scores, i.e.
        :math:`S^{mz}` and :math:`S^{intensity}`. In both cases, each peak :math:`k` is scored,
        comparing the measured value :math:`i` with the calculated value
        :math:`j` (equation 1), whereas a perfect match is 1. Each peak of the
        isotopologue that has a relative intensity (relative to the maximum
        intensity isotope peak) :math:`r_{k}` above the matching threshold
        (by default 1% of the maximum intensity isotope peak) is matched and
        scored.

        .. math::
           :nowrap:

            \\begin{equation}
                s^{}_{ijk} \\in [0, 1]
            \\end{equation}


        **The m/z score**

        For each peak :math:`k`, the m/z similarity between measured value
        :math:`i` and the calculated value :math:`j` is defined as

        .. math::
           :nowrap:

            \\begin{equation}
                s^{mz}_{ijk} = 1 - (\\frac{\\delta^{mz}_{ijk}}{\\alpha})
            \\end{equation}

        Whereas :math:`delta^{mz}_{ijk}` the difference in ppm between measured
        :math:`mz_{ik}` and calculated :math:`mz_{jk}` and :math:`\\alpha`
        defines the range in ppm, in which the score decreases from 1 to 0 in a
        linear fashion. In principle, :math:`\\alpha` is equal to the precision
        of the measurement defined by the user (pyQms parameter “REL_MZ_RANGE”,
        default 5 ppm, http://pyqms.readthedocs.io/en/latest/params.html).
        For example, if the difference between measured and theoretical m/z
        values would be 2.5 ppm, then the :math:`s^{mz}_{ijk}` score for this
        peak :math:`k` would be 0.5.

        The total m/z score for all peaks termed :math:`S^{mz}` is the weighted
        sum of all single similarity m/z scores :math:`s^{mz}_{ijk}` (equation 3).
        The weighting is defined by the theoretical intensity of the peak
        :math:`k` relative to the highest peak in the theoretical isotope
        pattern, termed :math:`r_{k}`.

        .. math::
           :nowrap:

            \\begin{equation}
                 S^{mz} = \\frac{\\sum\\limits_{}^k s^{mz}_{ijk} r_{k} }{\\sum\\limits_{}^k r_{k}}
            \\end{equation}


        **The intensity score**

        Prior to intensity scoring, the scaling factor :math:`\\sigma` is
        calculated by comparing the intensities of the measured :math:`i` and
        calculated :math:`j` intensities for all peaks :math:`k` within the
        matching threshold (see above). This scaling factor is calculated by
        dividing the weighted sum of the measured intensity by the weighted sum
        of the theoretical intensities (equation 4).

        .. math::
           :nowrap:

            \\begin{equation}
                 \\sigma = \\frac{\\sum\\limits_{}^k intensity_{ik} r_{k} }{\\sum\\limits_{}^k intensity_{jk} r_{k}}
            \\end{equation}

        Using this scaling factor, which is equal to the abundance of the
        measured molecule, one can calculate :math:`\\delta^{intensity}_{ijk}`,
        which is the relative intensity error between measured and theoretical
        intensity for each peak :math:`k` (equation 5).

        .. math::
           :nowrap:

            \\begin{equation}
                 \\delta^{intensity}_{ijk} = \\frac{ \\left|intensity_{ik} - \\sigma intensity_{jk}\\right|}{\\sigma intensity_{jk}}
            \\end{equation}

        The intensity score of peak :math:`k` is then defined (equation 6).

        .. math::
           :nowrap:

            \\begin{equation}
                 s^{intensity}_{ijk} = 1 - (\\frac{\\delta^{intensity}_{ijk}}{1 - r_{k} + \\epsilon })
            \\end{equation}

        In analogy to the m/z score (:math:`s^{mz}_{ijk}`), the denominator
        defines the range in which the peak based intensity score decreases
        from 1 to 0. However, in contrast to the m/z score, the intensity error
        has to be weighted by the abundance of each peak (1 - :math:`r_{k}` )
        as more abundant peaks can be measured more accurately than smaller
        peaks. Additionally, we introduced ϵ (pyQms parameter “REL_I_RANGE”,
        default 0.2), which represents the most conservative relative error
        applied to the most precisely measured peak (:math:`r_{k}` = 1). Thus,
        the overall relative error (denominator) will increase with lower peaks
        The total intensity score :math:`S^{intensity}` is the weighted sum of
        all similarity scores :math:`k` in analogy to the :math:`S^{mz}` score:


        .. math::
           :nowrap:

            \\begin{equation}
                S^{intensity} = \\frac{\\sum\\limits_{}^k s^{intensity}_{ijk} r_{k} }{\\sum\\limits_{}^k r_{k}}
            \\end{equation}

        **The combined final score: mScore**
        The final score is termed mScore and is a sum of :math:`S^{mz}` and
        :math:`S^{intensity}`. However, because some machines can measure m/z
        much more accurately then intensities, we introduced :math:`\\xi` to
        allow for flexibilities depending on the type of mass spectrometer used.
        :math:`\\xi` (the pyQms parameter “MZ_SCORE_PERCENTILE”, default 0.4)
        is the fraction the :math:`S^{mz}` score is weighted into the sum.
        Thus, the final mScore is defined as:


        .. math::
           :nowrap:

            \\begin{equation}
                mScore = \\xi S^{mz} + (1 - \\xi) S^{intensity}
            \\end{equation}



        ..
           commented out, this is the compacted form

           S^{mz} = \\frac{
                    \\sum\\limits_{}^k
                        \\big(1 - (
                            \\frac{
                                \\left| m^{mz}_{k} - c^{mz}_{k} \\right|
                            }
                            {
                                c^{mz}_{k} \\alpha
                            }
                        )\\big) \\delta(k)
                }
                {
                    \\sum\\limits_{}^k \\delta(k)
            }
            \\\\

            \\xi = \\frac{
                    \\sum\\limits_{}^k m^{intensity}_{k} \\delta(k)
                }
                {
                    \\sum\\limits_{}^k c^{intensity}_{k} \\delta(k)
                }
            \\\\\

            S^{intensity} = \\frac{
                \\sum\\limits_{}^k
                    \\big(1 - (
                        \\frac{
                            \\left| m^{intensity}_{k} - \\xi c^{intensity}_{k} \\right|
                        }
                        {
                            \\xi c^{intensity}_{k} (
                            1 - \\delta(k) + \\omega )
                        }
                    )\\big) \\delta(k)
            }
            {
                \\sum\\limits_{}^k \\delta(k)
            }
            \\\\\
            S = \\alpha
                \\bigg(
                         \\frac{
                                \\sum\\limits_{}^k
                                    \\big(1 - (
                                        \\frac{
                                            \\left| m^{mz}_{k} - c^{mz}_{k} \\right|
                                        }
                                        {
                                            c^{mz}_{k} \\alpha
                                        }
                                    )\\big) \\delta(k)
                            }
                            {
                                \\sum\\limits_{}^k \\delta(k)
                            }
                \\bigg)
                + ( 1 - \\alpha )
                \\bigg(
                        \\frac{
                                \\sum\\limits_{}^k
                                    \\big(1 - (
                                        \\frac{
                                            \\left| m^{intensity}_{k} - \\xi c^{intensity}_{k} \\right|
                                        }
                                        {
                                            \\xi c^{intensity}_{k} (
                                            1 - \\delta(k) + \\omega )
                                        }
                                    )\\big) \\delta(k)
                            }
                            {
                                \\sum\\limits_{}^k \\delta(k)
                            }
                \\bigg)


        """
        # print('-'*40)
        # if min_abundance is None:
        #     min_abundance = self.params['MOLECULE_MIN_ABUNDANCE']
        total_measured_i = 0
        total_calculated_i = 0
        ri_sum = 0
        for mmz, mi, ri, cmz, ci in matched_peaks:
            if ri < self.params["MIN_REL_PEAK_INTENSITY_FOR_MATCHING"]:
                continue
            if mmz is not None:
                total_calculated_i += ci * ri
                total_measured_i += mi * ri
            ri_sum += ri  # math.log(1+ri,10)

        scaling_factor = total_measured_i / float(total_calculated_i)
        # NOTE: old_pyQms scaling is not weighted by rel_i
        mz_score = 0
        i_score = 0
        mz_range = self.params["REL_MZ_RANGE"]
        i_range = self.params["REL_I_RANGE"]
        for mmz, mi, ri, cmz, ci in matched_peaks:
            if ri < self.params["MIN_REL_PEAK_INTENSITY_FOR_MATCHING"]:
                continue
            if mmz is None:
                continue
            si = ci * scaling_factor
            if cmz > sys.float_info.epsilon:
                rel_mz_error = abs(mmz - cmz) / cmz
                if rel_mz_error > mz_range:
                    local_mz_score = 0
                else:
                    local_mz_score = 1 - (rel_mz_error / mz_range)
                mz_score += local_mz_score * ri  # math.log(1+ri,10)#* ri

            if si > sys.float_info.epsilon:
                rel_i_error = abs(mi - si) / si
                scaled_i_range = 1.0 + i_range - ri
                if rel_i_error > scaled_i_range:
                    local_i_score = 0
                else:
                    local_i_score = 1 - (rel_i_error / scaled_i_range)
                i_score += local_i_score * ri  # math.log(1+ri,10) #* ri
        #
        mz_score = mz_score / ri_sum
        mz_score_fraction = mz_score * mz_score_percentile
        i_score = i_score / ri_sum
        i_score_fraction = i_score * (1 - mz_score_percentile)
        score = mz_score_fraction + i_score_fraction
        # if score > 1:
        #     print(matched_peaks,'WW')
        #     print( scaling_factor )
        #     print( i_score, mz_score)
        #     print([ri for mmz,mi,ri,cmz,ci in matched_peaks])
        #     print(sum([ri for mmz,mi,ri,cmz,ci in matched_peaks]))
        #     print(tmp)
        #     print(sum(tmp))
        #     exit(1)
        # print('score:\t{0}'.format(score))
        # print('scaling:',scaling_factor)
        # print()
        return score, scaling_factor

    def _slice_list(self, source_list, borders, tolerance=2):
        """
        Bisect a given list by using the minimum and maximum of a defined border
        (list or tuple of values)
        """
        lower_value = list(borders[0])
        upper_value = list(borders[1])
        is_numpy_array = getattr(source_list, "tolist", False)
        if is_numpy_array is False:
            # normal arrays ...
            try:
                min_pos = bisect.bisect(source_list, lower_value) - tolerance
                max_pos = bisect.bisect(source_list, upper_value) + tolerance
            except:
                # the list has tuples instead of list
                min_pos = bisect.bisect(source_list, tuple(lower_value)) - tolerance
                max_pos = bisect.bisect(source_list, tuple(upper_value)) + tolerance

            if min_pos < 0:
                min_pos = 0
            if max_pos > len(source_list):
                max_pos = len(source_list)

            r_list = source_list[min_pos:max_pos]
        else:
            r_list = source_list[
                np.where(
                    (lower_value[0] - tolerance < source_list[:, 0])
                    * (source_list[:, 0] < upper_value[0] + tolerance)
                )
            ]
        return r_list

    def _transform_mz_to_set(self, mz):
        """
        Internal function which generates a set of transformed mz values.

        An error is included based on the measured precision.
        The measured precision (e.g. 5ppm or 5e-6 Da ) is used to define lower
        and upper borders.

        All mz values are transformed using an internal factor
        (e.g. 1000 ) to use integers instead of floats. This results in a faster
        processing in the algorithm.
        """
        mz_error = mz * self.params["REL_MZ_RANGE"]
        lower_mz = mz - mz_error
        tlower_mz = lower_mz * self.params["INTERNAL_PRECISION"]
        upper_mz = mz + mz_error
        tupper_mz = upper_mz * self.params["INTERNAL_PRECISION"]

        tmp = set(range(int(round(tlower_mz)), int(round(tupper_mz)) + 1))
        return tmp

    def _transform_spectrum(self, mz_i_list, mz_range=None):
        """
        Internal function to transform a spectrum.

        Returns a set of transformed mz values from a given
        spectrum (e.g. mz,intensity as tuples in a list). The mz values are
        transformed by the INTERNAL_PRECISION (i.e. 1000) to use integers
        instead of floats.
        Additionally a lookup is returned which maps the transformed values to
        the original m/z values including their intensities.
        """
        tmz_set = set()
        tmz_lookup = {}
        if mz_range is None:
            target_mz_list = mz_i_list
        else:
            target_mz_list = self._slice_list(mz_i_list, [(x, 0.001) for x in mz_range])
        # for mz, intensity in target_mz_list:
        #     tmz = int(round( mz * self.params['INTERNAL_PRECISION'] ))
        #     tmz_set.add( tmz )
        #     try:
        #         tmz_lookup[ tmz ].append( (mz, intensity) )
        #     except:
        #         tmz_lookup[ tmz ] = [ (mz, intensity) ]

        is_numpy_array = getattr(target_mz_list, "tolist", False)
        if is_numpy_array is False:
            for mz, intensity in target_mz_list:
                tmz = int(round(mz * self.params["INTERNAL_PRECISION"]))
                tmz_set.add(tmz)
                try:
                    tmz_lookup[tmz].append((mz, intensity))
                except:
                    tmz_lookup[tmz] = [(mz, intensity)]
        else:
            tmz = np.round(
                target_mz_list[:, 0] * self.params["INTERNAL_PRECISION"]
            ).astype(int)
            for pos, tmz_entry in enumerate(tmz):
                try:
                    tmz_lookup[tmz_entry].append(
                        (target_mz_list[pos][0], target_mz_list[pos][1])
                    )
                except:
                    tmz_lookup[tmz_entry] = [
                        (target_mz_list[pos][0], target_mz_list[pos][1])
                    ]
            tmz_set = set(tmz)

        return tmz_set, tmz_lookup


if __name__ == "__main__":
    print(__doc__)
