"""
PyHELM - A class for working with HELM

Copyright 2014, quattro research GmbH
"""

import re


class NotationError (Exception):
    """Custom exception for invalid HELM"""


# pylint: disable=redefined-builtin
# noinspection PyShadowingBuiltins
class Polymer:
    """Class representing a polymer chain"""
    def __init__(self, type, data='', name=''):
        self.type = type
        self.data = data
        self.name = name


# noinspection PyShadowingNames
class Connection:
    """  A connection has 3 Attributes: Polymer, Start, Stop; All of them are tuples with two elements

    Example: (Polymer1,Polymer2),(20,'R2'),(3,'R1') or (Polymer1,Polymer2),(10,'pair'),(3,'pair')
    """
    def __init__(self, poly, start, stop):
        self.polymers = poly
        self.start = start
        self.stop = stop

    def validate(self):
        """Validate connection between polymer chains"""
        if len(self.polymers) == 2 and len(self.start) == 2 and len(self.stop) == 2:
            if self.polymers[0].countmonomer() > self.start[0] and self.polymers[1].countmonomer() > self.stop[0]:
                return True
        return False


# noinspection PyUnresolvedReferences,RegExpRedundantEscape
class HelmObj:
    """ A HelmObj is a complete polymer made of simple polymers, connections and attributes, separated by the '$' sign.

    self.polymers = List of polymer objects
    self.connections = List of connections
    self.attributes = List of (Polymer,String)
    self.notation = HELM String
    """

    def __init__(self):
        self.polymers = []
        self.connections = []
        self.attributes = []
        self.data = ''

# ####        self.monomerDB = os.path.dirname(__file__) + '/MonomerDB.xml'
# ####
# ####        if not os.path.isfile(self.monomerDB):
# ####            print('could not find monomerDB')
# ####            return
# ####
# ####        self.tree = ET.parse(self.monomerDB)
# ####        polymers = self.tree.findall("PolymerList/Polymer")
#        for p in polymers:
#            print(p.attrib['polymerType'])
#
#        print('monomerDB initialized')

    # noinspection PyPep8Naming,PyMethodMayBeStatic
    def validate_helm(self, helm):  # pylint: disable=no-self-use
        """Validate input HELM string"""
        # #check $ signs
        sections = re.split(r'\$', helm)
        if len(sections) != 5:
            return False

        # #check polymer connections
        first_list = re.findall(r'\w+(?=\{)', sections[0])  # all polymers from first element
        second_list = re.findall(r'\w+(?=[\{,])', sections[1])   # all polymers from second element
        for element in second_list:
            if element not in first_list:
                return False
        return True

    # noinspection PyPep8Naming,PyShadowingNames
    def parse_helm(self, helm):
        """Parse input HELM string into objects of Polymers, Connection and Attributes"""
        # clear the list of objects
        del self.polymers[:]
        del self.connections[:]
        del self.attributes[:]

        if not self.validate_helm(helm):
            raise NotationError('Invalid HELM notation')

        sections = re.split(r'\$', helm)
        poly_section = sections[0]
        poly_strings = re.split(r'\|', poly_section)
        for p in poly_strings:
            p_name = re.match(r'\w+(?=\{)', p).group()
            p_type = re.match(r'\D+', p_name).group()
            p_data = p[len(p_name) + 1:len(p) - 1]
            self.polymers.append(Polymer(p_type, p_data, p_name))

        if sections[1] != '':
            conn_section = sections[1]
            conn_strings = re.split(r'\|', conn_section)
            for c in conn_strings:
                items = re.split(r'[,:-]', c)
                polys = [poly for poly in self.polymers if poly.name == items[0]]
                polys.extend([poly for poly in self.polymers if poly.name == items[1]])
                start = (int(items[2]), items[3])
                stop = (int(items[4]), items[5])
                self.connections.append(Connection(polys, start, stop))

    # noinspection PyPep8Naming,PyMethodMayBeStatic
    def get_polymer_names_from_helm(self, helm):  # pylint: disable=no-self-use
        """Get Polymer chain names from a HELM string"""
        sections = re.split(r'\$', helm)
        return re.findall(r'\w+(?=\{)', sections[0])

    # noinspection PyPep8Naming,PyMethodMayBeStatic
    def get_polymers_from_helm(self, helm):  # pylint: disable=no-self-use
        """Get Polymer chain name and sequence from a HELM string"""
        sections = re.split(r'\$', helm)
        poly_section = sections[0]
        return re.split(r'\|', poly_section)

    # noinspection PyPep8Naming,PyMethodMayBeStatic
    def get_monomers_from_polymer(self, helm, polymer_name):  # pylint: disable=no-self-use,inconsistent-return-statements
        """Get all monomers from a Polymer chain - the sequence"""
        sections = re.split(r'\$', helm)
        poly_section = sections[0]
        polymers = re.split(r'\|', poly_section)
        for p in polymers:
            if p.find(polymer_name) >= 0:
                # remove name and curly brackets
                mono_str = p[len(polymer_name) + 1:len(p) - 1]
                return re.split(r'\.', mono_str)

    # noinspection PyPep8Naming
    def get_distinct_monomers_from_polymer(self, helm, polymer_name):
        """Return a set of unique monomers in a polymer chain"""
        all_monomers = self.get_monomers_from_polymer(helm, polymer_name)
        return set(all_monomers)


if __name__ == '__main__':
    # example use
    my_helm = HelmObj()

    my_helm.parse_helm('RNA1{[am6]P.R(C)P.R(U)P.R(U)P.R(G)P.R(A)P.R(G)P.R(G)}|PEPTIDE1{A.C.G.K.E.D.K.R}|CHEM1{SMCC}$PEPTIDE1,CHEM1,2:R3-1:R2|RNA1,CHEM1,1:R1-1:R1$$$')  # pylint: disable=line-too-long
    for polymer in my_helm.polymers:
        print(polymer.type)
        print(polymer.data)
        print('')

    for conn in my_helm.connections:
        print(conn.polymers[0].name + str(conn.start) + ' -> ' + conn.polymers[1].name + str(conn.stop))
