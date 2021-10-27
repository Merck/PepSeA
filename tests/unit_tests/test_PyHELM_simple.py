import pytest

from alignment.PyHELM_simple import HelmObj, NotationError, Polymer, Connection


# Pytest fixtures
@pytest.fixture
def helm_obj():
    # Return empty HELM object
    return HelmObj()


@pytest.fixture
def helm_string():
    helm_sequence = "PEPTIDE1{[Ac].[dApe].P.L.Y.I.S.Y.D.P.V.[Ape].[NH2]}|PEPTIDE2{[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]}|CHEM1{[PEG2diacid]}$PEPTIDE1,CHEM1,18:R3-1:R2|PEPTIDE2,CHEM1,18:R3-1:R1|PEPTIDE2,PEPTIDE2,1:R1-15:R3|PEPTIDE1,PEPTIDE1,1:R1-15:R3$$$V2.0"

    return helm_sequence


"""
Test the correct raise of NotationError
"""


def test_notation_error():
    with pytest.raises(Exception) as exception:
        raise NotationError()

    assert exception.typename == "NotationError"


"""
Test correct initiation of the Polymer class
"""


def test_polymer_init():

    polymer = Polymer("PEPTIDE", "[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]", "PEPTIDE1")

    assert polymer.type == "PEPTIDE"
    assert polymer.data == "[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]"
    assert polymer.name == "PEPTIDE1"


def test_polymer_init_incorrectly():
    """ Test the initiation with no arguments"""

    with pytest.raises(Exception) as exception:
        Polymer()

    assert exception.typename == "TypeError"


"""
Test correct initiation of the Connection class
"""


def test_connection_init():
    polymers = ("PEPETIDE1", "CHEM1")
    start = (18, "R3")
    finish = (1, "R1")

    connection = Connection(polymers, start, finish)

    assert connection.polymers == polymers
    assert connection.start == start
    assert connection.stop == finish


def test_connection_init_incorrectly():
    """ Test the initiation with no arguments"""

    with pytest.raises(Exception) as exception:
        Connection()

    assert exception.typename == "TypeError"


"""
Below are unit tests of the methods of HelmObj class
"""


def test_helmobj_init():
    """Test the correct initiation"""

    helm = HelmObj()

    assert helm.polymers == []
    assert helm.connections == []
    assert helm.attributes == []
    assert helm.data == ""


"""
Below are tests for validate_helm() function
"""


def test_validate_helm(helm_obj, helm_string):
    """Test validation of the correct HELM string"""

    result = helm_obj.validate_helm(helm_string)
    assert result is True


def test_validate_helm_wrong_number_of_sections(helm_obj):
    """Test the validation of incorrect HELM string - incorrect number of sections (parts of HELM, separated by $ symbol)"""

    helm_string = 'PEPTIDE1{[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]}|PEPTIDE2{[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]}|CHEM1{[PEG2diacid]}$PEPTIDE1,CHEM1,18:R3-1:R2|PEPTIDE2,CHEM1,18:R3-1:R1|PEPTIDE2,PEPTIDE2,1:R1-15:R3|PEPTIDE1,PEPTIDE1,1:R1-15:R3$$'
    result = helm_obj.validate_helm(helm_string)

    assert result is False


# Fix function!
def test_validate_helm_wrong_connections(helm_obj):

    helm_string = 'PEPTIDE1{[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]}|PEPTIDE2{[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]}|CHEM1{[PEG2diacid]}$PEPTIDE1,CHEM1,18:R3-1:R2|PEPTIDE2,CHEM1,18:R3-1:R1|PEPTIDE2,PEPTIDE3,1:R1-15:R3|PEPTIDE1,PEPTIDE1,1:R1-15:R3$$$V2.0'
    validation = helm_obj.validate_helm(helm_string)
    assert validation is False


"""
Below are tests for parse_helm() function
"""


def test_parse_helm(helm_obj, helm_string):

    # Parse the string
    helm_obj.parse_helm(helm_string)

    # Manually define the polymers, which are present in the input HELM string
    polymer1 = Polymer("PEPTIDE", "[Ac].[dApe].P.L.Y.I.S.Y.D.P.V.[Ape].[NH2]", "PEPTIDE1")
    polymer2 = Polymer("PEPTIDE", "[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]", "PEPTIDE2")
    polymer3 = Polymer("CHEM", "[PEG2diacid]", "CHEM1")

    # List of polymers
    polymers = [polymer1, polymer2, polymer3]

    # Manually define the connections between polymers, which are present in the input HELM string
    connection1 = Connection([polymer1, polymer3], (18, "R3"), (1, "R2"))
    connection2 = Connection([polymer2, polymer3], (18, "R3"), (1, "R1"))
    connection3 = Connection([polymer2, polymer2], (1, "R1"), (15, "R3"))
    connection4 = Connection([polymer1, polymer1], (1, "R1"), (15, "R3"))

    # List of connections
    connections = [connection1, connection2, connection3, connection4]

    # Iterate through parsed and manually defined polymers and compare them
    for poly_parsed, poly_expected in zip(helm_obj.polymers, polymers):
        # Compare the name, type and data of polymers

        assert poly_parsed.type == poly_expected.type
        assert poly_parsed.name == poly_expected.name
        assert poly_parsed.data == poly_expected.data

    # Iterate through the parsed and manually defined connections and compare them
    for i, conn in enumerate(zip(helm_obj.connections, connections)):

        con_parsed, con_expected = conn

        # Iterate through the polymers of parsed connections, and manually defined connections (which are accessed by index)
        for poly_parsed, poly_expected in zip(con_parsed.polymers, connections[i].polymers):
            # Compare the name, type and data of polymers

            assert poly_parsed.name == poly_expected.name
            assert poly_parsed.type == poly_expected.type
            assert poly_parsed.data == poly_expected.data

        # Compare start and stop variables of the connection
        assert con_parsed.start == con_expected.start
        assert con_parsed.stop == con_expected.stop


def test_parse_helm_invalid_helm(helm_obj):
    """Test that exception is raised, given the corrupted HELM string"""

    # Corrupted HELM string
    helm_string = 'PEPTIDE1{[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]}|PEPTIDE2{[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]}|CHEM1{[PEG2diacid]}$PEPTIDE1,CHEM1,18:R3-1:R2|PEPTIDE2,CHEM1,18:R3-1:R1|PEPTIDE2,PEPTIDE2,1:R1-15:R3|PEPTIDE1,PEPTIDE1,1:R1-15:R3$$'

    with pytest.raises(Exception) as exception:
        helm_obj.parse_helm(helm_string)

    assert str(exception.value) == "Invalid HELM notation"


"""
Below are tests for get_polymer_names_from_helm() method
"""


def test_get_polymer_names_from_helm(helm_obj, helm_string):
    # Test extraction of names of peptides presented in HELM sequence

    result = helm_obj.get_polymer_names_from_helm(helm_string)
    expected = ["PEPTIDE1", "PEPTIDE2", "CHEM1"]
    assert result == expected


"""
Below are tests for get_polymers_from_helm() method
"""


def test_get_polymers_from_helm(helm_obj, helm_string):
    # Test extraction of distinct polymers from HELM string

    polymer1 = "PEPTIDE1{[Ac].[dApe].P.L.Y.I.S.Y.D.P.V.[Ape].[NH2]}"
    polymer2 = "PEPTIDE2{[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]}"
    polymer3 = "CHEM1{[PEG2diacid]}"

    expected_polymers = [polymer1, polymer2, polymer3]

    polymers = helm_obj.get_polymers_from_helm(helm_string)

    assert polymers == expected_polymers


"""
Below are tests for get_polymers_from_helm() method
"""


def test_get_monomers_from_polymer(helm_obj, helm_string):
    # Test the extraction of monomers for each distinct peptide, presented in HELM sequence

    data = {"PEPTIDE1": ["[Ac]", "[dApe]", "P", "L", "Y", "I", "S", "Y", "D", "P", "V", "[Ape]", "[NH2]"],
            "PEPTIDE2": ["[ClAc]", "F", "S", "V", "[Sar]", "[Ahp]", "R", "R", "[NMeF]", "V", "A", "[Ahp]", "S", "[Bip]", "C", "G", "G", "K", "[NH2]"],
            "CHEM1": ["[PEG2diacid]"]}

    for key in data:
        monos = helm_obj.get_monomers_from_polymer(helm_string, key)
        assert monos == data[key]


"""
Below are tests for get_distinct_monomers_from_polymer() method
"""


def test_get_distinct_monomers_from_polymer(helm_obj, helm_string):
    monos = helm_obj.get_distinct_monomers_from_polymer(helm_string, "PEPTIDE1")
    expected_monos = {'V', '[NH2]', '[dApe]', 'L', 'Y', 'D', '[Ape]', 'S', '[Ac]', 'I', 'P'}

    assert monos == expected_monos
