import json
import pytest

from alignment.ApiUtils import extract_helm_from_json, extract_subpeptide, extract_monomer, \
    extract_aligned_sequences, json_output, escape_double_quotes_in_input


# Pytest fixtures
@pytest.fixture
def helm_seq():
    """Returns HELM string"""
    return "PEPTIDE1{[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]}|PEPTIDE2{[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]}|CHEM1{[PEG2diacid]}$PEPTIDE1,CHEM1,18:R3-1:R2|PEPTIDE2,CHEM1,18:R3-1:R1|PEPTIDE2,PEPTIDE2,1:R1-15:R3|PEPTIDE1,PEPTIDE1,1:R1-15:R3$$$V2.0"


@pytest.fixture
def path_to_data():
    # Returns path to data

    path = "tests/unit_tests/resources/test_ApiUtils/"
    return path


"""
Below are tests for function extract_helm_from_json()
"""


def test_extract_helm_from_json(path_to_data):
    """Test extraction of HELM strings from JSON input"""

    path = path_to_data + 'example.json'
    with open(path) as file:
        input_json = json.load(file)
        output = extract_helm_from_json(input_json)
        expected_output = "PEPTIDE1{F.C.R.C.A.C.C.D.M.K.L}|PEPTIDE2{F.C.R.C.A.C.C.D.M.K.a" \
                          "MeL}$$$$\nPEPTIDE1{F.C.R.C.A.C.C.D.M.K.L}|PEPTIDE2{F.C.R.C.A.C.C.D.M.K.aMeL}$$$$"

        assert output == expected_output


@pytest.mark.xfail(raises=KeyError)
def test_extract_helm_from_json_fail():
    """Test that KeyError is raised"""

    input_helm = [{"ID": "L-000000000"},
                  {"ID": "L-000000001"}]

    # Call function on input with no HELM entries
    extract_helm_from_json(input_helm)


"""
Below are tests for function extract_subpeptide()
"""


def test_extract_subpeptide(helm_seq):
    """Test the extraction of subpeptide from HELM string"""

    expected_output = "[ClAc].F.S.V.[Sar].[Ahp].R.R.[NMeF].V.A.[Ahp].S.[Bip].C.G.G.K.[NH2]"
    output = extract_subpeptide("PEPTIDE1", helm_seq)

    assert output == expected_output


def test_extract_subpeptide_missing(helm_seq):
    """Test that function returns None, whet peptide is not in the HEL sequence"""

    output = extract_subpeptide("PEPTIDE3", helm_seq)

    assert output is None


"""
Below are tests for function extract_monomer()
"""


def test_extract_monomer():
    """Test the correct work of generator function extract_monomer()"""

    aligned_seq = "[Ac]--------------[LysN3]RRRRCPLYIS[NMeY]DPVCRRRR[NH2]"
    monos = ['Ac', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'LysN3', 'R', 'R', 'R', 'R', 'C', 'P', 'L', 'Y', 'I', 'S', 'NMeY', 'D', 'P', 'V', 'C', 'R', 'R', 'R', 'R', 'NH2']

    for extracted, mono in zip(extract_monomer(aligned_seq), monos):
        assert extracted == mono


def test_extract_monomer_corrupted_seq():
    """Test the function on corrupted sequence raises assertion error"""

    # In the following sequence, nnAA LysN3 is missing opening bracket. Should it lead to error?
    aligned_seq = "[Ac]----[LysN3RRCPLYISNMeYDPRNH2"

    with pytest.raises(Exception) as excinfo:
        [el for el in extract_monomer(aligned_seq)]

    expected_message = f"Non-valid notation of non-natural amino acids in the following sequence {aligned_seq}"
    assert str(excinfo.value) == expected_message


"""
Below are tests for function extract_aligned_sequences()
"""


def test_extract_aligned_sequences(path_to_data):

    path = path_to_data + "aligned.json"
    with open(path) as file:
        input_json = json.load(file)
        expected_output = "> PEPTIDE2\n[ClAc]-FYSW[Sar]NYWSYY[Sar]WC-[NH2]\n> PEPTIDE2\n[ClAc]-FYSW[Sar]NYWSYA[Sar]WCG[NH2]\n> PEPTIDE2\n[ClAc]-FYSW[Sar]NYWSYYAWCG[NH2]\n> PEPTIDE2\n[ClAc]-FYSW[Sar]NYWSYY[Sar]ACG[NH2]"
        output = extract_aligned_sequences(input_json)

        assert output == expected_output


"""
Below are tests for json_output() method
"""


def test_json_output(path_to_data):
    """Test the correct creation of JSON output"""

    aligned_path = path_to_data + "raw_output.json"
    input_path = path_to_data + "raw_input.json"
    expected_path = path_to_data + "json_output.json"

    with open(aligned_path) as aligned, open(input_path) as input_json, open(expected_path) as out:
        aligned_data = json.load(aligned)
        input_data = json.load(input_json)
        expected_output = json.load(out)

        output = json_output(aligned_data, input_data)

        assert output == expected_output


"""
Below are test for escape_double_quotes_in_input() method
"""


def test_escape_double_quotes_in_input():
    """Test the correct addition of escape characters to the HELM string"""

    input_data = [{'ID': 'L-006579559-000J001',
                   'HELM': 'PEPTIDE1{[ClAc].[dF].[dV].[dApe].[dW].[dY].[dD].[dNMeF].G.[dL].[dApe].[dL].[dC].[NH2]}|CHEM1{[tRCM]}$PEPTIDE1,CHEM1,4:R3-1:R1|PEPTIDE1,CHEM1,11:R3-1:R2|PEPTIDE1,PEPTIDE1,1:R1-13:R3$${"chiral":"ena"}$V2.0'
                   },
                  {'ID': 'L-006579556-000H001',
                   'HELM': 'PEPTIDE1{[ClAc].[dF].[dV].[dY].[rbaFW].[dY].[dD].[dNMeF].G.[dApe].[dD].[dL].[dC].[NH2]}|CHEM1{[cRCM]}$PEPTIDE1,PEPTIDE1,1:R1-13:R3|PEPTIDE1,CHEM1,5:R3-1:R1|PEPTIDE1,CHEM1,10:R3-1:R2$${"chiral":"mix"}$V2.0'
                   }]

    expected_output = [{'ID': 'L-006579559-000J001',
                        'HELM': 'PEPTIDE1{[ClAc].[dF].[dV].[dApe].[dW].[dY].[dD].[dNMeF].G.[dL].[dApe].[dL].[dC].[NH2]}|CHEM1{[tRCM]}$PEPTIDE1,CHEM1,4:R3-1:R1|PEPTIDE1,CHEM1,11:R3-1:R2|PEPTIDE1,PEPTIDE1,1:R1-13:R3$${\\\"chiral\\\":\\\"ena\\\"}$V2.0'
                        },
                       {'ID': 'L-006579556-000H001',
                       'HELM': 'PEPTIDE1{[ClAc].[dF].[dV].[dY].[rbaFW].[dY].[dD].[dNMeF].G.[dApe].[dD].[dL].[dC].[NH2]}|CHEM1{[cRCM]}$PEPTIDE1,PEPTIDE1,1:R1-13:R3|PEPTIDE1,CHEM1,5:R3-1:R1|PEPTIDE1,CHEM1,10:R3-1:R2$${\\\"chiral\\\":\\\"mix\\\"}$V2.0'
                        }]

    result = escape_double_quotes_in_input(input_data)

    for output, expected in zip(result, expected_output):

        assert output["HELM"] == expected["HELM"]
