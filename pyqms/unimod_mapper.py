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
import os
import codecs
import xml.etree.ElementTree as ET


class UnimodMapper(object):
    """
    UnimodMapper class that creates lookup to the unimod.xml located in
    kb/ext/unimod.xml and offers several helper methods.

    Mapping from e.g. name to composition or unimod ID to mass is possible.
    
    Please refer to `unimod`_ for further informations on modifications
    including naming, formulas, masses etc.

    .. _unimod:
        http://www.unimod.org/modifications_list.php?
    
    """

    def __init__(self):
        self._data_list = None
        self._mapper = None

        self.unimod_xml_name = "unimod.xml"
        # self.data_list = self._parseXML()
        # self.mapper    = self._initialize_mapper()

    @property
    def data_list(self):
        if self._data_list is None:
            self._data_list = self._parseXML()
        return self._data_list

    @data_list.setter
    def data_list(self, data_list):
        self._data_list = data_list
        return

    @property
    def mapper(self):
        if self._mapper is None:
            self._mapper = self._initialize_mapper()
        return self._mapper

    @mapper.setter
    def mapper(self, mapper):
        self._mapper = mapper
        return

    def _parseXML(self):
        is_frozen = getattr(sys, "frozen", False)
        if is_frozen:
            xmlFile = os.path.normpath(
                os.path.join(os.path.dirname(sys.executable), self.unimod_xml_name)
            )
        else:
            xmlFile = os.path.normpath(
                os.path.join(
                    os.path.dirname(__file__), "kb", "ext", self.unimod_xml_name
                )
            )
        data_list = []
        print("> Parsing unimod.xml from {0}".format(xmlFile))
        if os.path.exists(xmlFile):
            unimodXML = ET.iterparse(
                codecs.open(xmlFile, "r", encoding="utf8"), events=(b"start", b"end")
            )
            collect_element = False
            for event, element in unimodXML:
                if event == b"start":
                    if element.tag.endswith("}mod"):
                        tmp = {
                            "unimodID": element.attrib["record_id"],
                            "unimodname": element.attrib["title"],
                            "element": {},
                            "specificity_sites": [],
                        }
                    elif element.tag.endswith("}delta"):
                        collect_element = True
                        tmp["mono_mass"] = float(element.attrib["mono_mass"])
                    elif element.tag.endswith("}element"):
                        if collect_element is True:
                            number = int(element.attrib["number"])
                            if number != 0:
                                tmp["element"][element.attrib["symbol"]] = number
                    elif element.tag.endswith("}specificity"):
                        amino_acid = element.attrib["site"]
                        if element.attrib["classification"] != "Artefact":
                            tmp["specificity_sites"].append(amino_acid)
                    else:
                        pass
                else:
                    # end element
                    if element.tag.endswith("}delta"):
                        collect_element = False
                    elif element.tag.endswith("}mod"):
                        data_list.append(tmp)
                    else:
                        pass

        else:
            print("No unimod.xml file found. Expected at {0}".format(xmlFile))
            sys.exit(1)
        return data_list

    def _initialize_mapper(self):
        """
        Set up the mapper, generates the index dict
        """
        mapper = {}
        for index, unimod_data_dict in enumerate(self.data_list):
            for key, value in unimod_data_dict.items():
                if key == "element":
                    MAJORS = ["C", "H"]
                    hill_notation = ""
                    for major in MAJORS:
                        if major in unimod_data_dict[key].keys():
                            hill_notation += "{0}({1})".format(
                                major, unimod_data_dict[key][major]
                            )
                    for symbol, number in sorted(unimod_data_dict[key].items()):
                        if symbol in MAJORS:
                            continue
                        hill_notation += "{0}({1})".format(symbol, number)

                    if hill_notation not in mapper.keys():
                        mapper[hill_notation] = []
                    mapper[hill_notation].append(index)
                elif key == "mono_mass":
                    if value not in mapper.keys():
                        mapper[value] = []
                    mapper[value].append(index)
                elif key == "specificity_sites":
                    pass
                else:
                    if value not in mapper.keys():
                        mapper[value] = index
        return mapper

    # name 2 ....
    def name2mass(self, unimod_name):
        """
        Converts unimod name to unimod mono isotopic mass

        Args:
            unimod_name (str): name of modification (as named in unimod)

        Returns:
            float: Unimod mono isotopic mass
        """
        return self._map_key_2_index_2_value(unimod_name, "mono_mass")

    def name2composition(self, unimod_name):
        """
        Converts unimod name to unimod composition

        Args:
            unimod_name (str): name of modification (as named in unimod)

        Returns:
            dict: Unimod elemental composition
        """
        return self._map_key_2_index_2_value(unimod_name, "element")

    def name2id(self, unimod_name):
        """
        Converts unimod name to unimod ID

        Args:
            unimod_name (str): name of modification (as named in unimod)

        Returns:
            int: Unimod ID
        """
        return self._map_key_2_index_2_value(unimod_name, "unimodID")

    def name2specificity_site_list(self, unimod_name):
        """
        Converts unimod name to list of specified amino acids or sites

        Args:
            unimod_name (str): name of modification (as named in unimod)

        Returns:
            list: list of specificity sites
        """
        list_2_return = None
        index = self.mapper.get(unimod_name, None)
        if index is not None:
            list_2_return = self._data_list_2_value(index, "specificity_sites")
        return list_2_return

    # unimodid 2 ....
    def id2mass(self, unimod_id):
        """
        Converts unimod ID to unimod mass

        Args:
            unimod_id (int): identifier of modification

        Returns:
            float: Unimod mono isotopic mass
        """
        return self._map_key_2_index_2_value(unimod_id, "mono_mass")

    def id2composition(self, unimod_id):
        """
        Converts unimod ID to unimod composition

        Args:
            unimod_id (int): identifier of modification

        Returns:
            dict: Unimod elemental composition
        """
        return self._map_key_2_index_2_value(unimod_id, "element")

    def id2name(self, unimod_id):
        """
        Converts unimod ID to unimod name

        Args:
            unimod_id (int): identifier of modification

        Returns:
            str: Unimod name
        """
        return self._map_key_2_index_2_value(unimod_id, "unimodname")

    # mass is ambigous therefore a list is returned
    def mass2name_list(self, mass):
        """
        Converts unimod mass to unimod name list,
        since a given mass can map to mutiple entries in the XML.

        Args:
            mass (float): mass of modification

        Returns:
            list: Unimod names
        """
        list_2_return = []
        for index in self.mapper[mass]:
            list_2_return.append(self._data_list_2_value(index, "unimodname"))
        return list_2_return

    def mass2id_list(self, mass):
        """
        Converts unimod mass to unimod name list,
        since a given mass can map to mutiple entries in the XML.

        Args:
            mass (float): mass of modification

        Returns:
            list: Unimod IDs
        """
        list_2_return = []
        index_list = self.mapper.get(mass, None)
        if index_list is not None:
            for index in index_list:
                list_2_return.append(self._data_list_2_value(index, "unimodID"))
        return list_2_return

    def mass2composition_list(self, mass):
        """
        Converts unimod mass to unimod element composition list,
        since a given mass can map to mutiple entries in the XML.

        Args:
            mass (float): mass of modification

        Returns:
            list: Unimod elemental compositions
        """

        list_2_return = []
        for index in self.mapper[mass]:
            list_2_return.append(self._data_list_2_value(index, "element"))
        return list_2_return

    def appMass2id_list(self, mass, decimal_places=2):
        """
        Creates a list of unimod IDs for a given approximate mass

        Args:
            mass (float): approximate mass of modification

        Keyword Arguments:
            decimal_places (int): Precision with which the masses in the
                Unimod is compared to the input, i.e. round( mass, decimal_places )

        Returns:
            list: Unimod IDs

        Examples::

            >>> import pyqms
            >>> U = pyqms.UnimodMapper()
            >>> U.appMass2id_list(18, decimal_places=0)
                ['127', '329', '608', '1079', '1167']

        """
        return_list = self._appMass2whatever(
            mass, decimal_places=decimal_places, entry_key="unimodID"
        )
        return return_list

    def appMass2element_list(self, mass, decimal_places=2):
        """
        Creates a list of element composition dicts for a given approximate mass

        Args:
            mass (float): approximate mass of modification

        Keyword Arguments:
            decimal_places (int): Precision with which the masses in the
                Unimod is compared to the input, i.e. round( mass, decimal_places )

        Returns:
            list: Dicts of elements

        Examples::

            >>> import pyqms
            >>> U = pyqms.UnimodMapper()
            >>> U.appMass2element_list(18, decimal_places=0)
                [{'F': 1, 'H': -1}, {'13C': 1, 'H': -1, '2H': 3},
                {'H': -2, 'C': -1, 'S': 1}, {'H': 2, 'C': 4, 'O': -2},
                {'H': -2, 'C': -1, 'O': 2}]


        """
        return_list = self._appMass2whatever(
            mass, decimal_places=decimal_places, entry_key="element"
        )
        return return_list

    def appMass2name_list(self, mass, decimal_places=2):
        """
        Creates a list of unimod names for a given approximate mass

        Args:
            mass (float): approximate mass of modification

        Keyword Arguments:
            decimal_places (int): Precision with which the masses in the
                Unimod is compared to the input, i.e. round( mass, decimal_places )

        Returns:
            list: Unimod names

        Examples::

            >>> import pyqms
            >>> U = pyqms.UnimodMapper()
            >>> U.appMass2name_list(18, decimal_places=0)
                ['Fluoro', 'Methyl:2H(3)13C(1)', 'Xle->Met', 'Glu->Phe', 'Pro->Asp']

        """
        return_list = self._appMass2whatever(
            mass, decimal_places=decimal_places, entry_key="unimodname"
        )
        return return_list

    def composition2name_list(self, composition):
        """
        Converts unimod composition to unimod name list,
        since a given composition can map to mutiple entries in the XML.

        Args:
            composition (dict): element composition (element, count pairs)

        Returns:
            list: Unimod names
        """
        list_2_return = []
        index_list = self.mapper.get(composition, None)
        if index_list is not None:
            for index in index_list:
                list_2_return.append(self._data_list_2_value(index, "unimodname"))
        return list_2_return

    def composition2id_list(self, composition):
        """
        Converts unimod composition to unimod name list,
        since a given composition can map to mutiple entries in the XML.

        Args:
            composition (dict): element composition (element, count pairs)

        Returns:
            list: Unimod IDs
        """
        list_2_return = []
        index_list = self.mapper.get(composition, None)
        if index_list is not None:
            for index in index_list:
                list_2_return.append(self._data_list_2_value(index, "unimodID"))
        return list_2_return

    def composition2mass(self, composition):
        """
        Converts unimod composition to unimod monoisotopic mass.

        Args:
            composition (dict): element composition (element, count pairs)

        Returns:
            float: monoisotopic mass
        """
        mass_2_return = None
        list_2_return = []
        index_list = self.mapper.get(composition, None)
        if index_list != None:
            for index in index_list:
                list_2_return.append(self._data_list_2_value(index, "mono_mass"))
            assert (
                len(set(list_2_return)) == 1
            ), """
            Unimod chemical composition {0}
            maps on different monoisotopic masses. This should not happen.
            """.format(
                composition
            )
            mass_2_return = list_2_return[0]
        return mass_2_return

    def _appMass2whatever(self, mass, decimal_places=2, entry_key=None):
        return_list = []
        for entry in self.data_list:
            umass = entry["mono_mass"]
            rounded_umass = round(float(umass), decimal_places)
            if abs(rounded_umass - mass) <= sys.float_info.epsilon:
                return_list.append(entry[entry_key])
        return return_list

    def _map_key_2_index_2_value(self, map_key, return_key):
        index = self.mapper.get(map_key, None)
        if index is None:
            print("Cant map", map_key, file=sys.stderr)
            return_value = None
        else:
            return_value = self._data_list_2_value(index, return_key)
        return return_value

    def _data_list_2_value(self, index, return_key):
        return self.data_list[index][return_key]


if __name__ == "__main__":
    print(__doc__)
    print(UnimodMapper.__doc__)
