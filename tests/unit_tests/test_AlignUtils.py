import os
import re
import json
import pytest

from alignment.AlignUtils import AlignUtils


# Pytest fixtures
@pytest.fixture
def path_to_test_data():
    """Returns path to the folder with test data"""
    return "tests/unit_tests/resources/test_AlignUtils/"


@pytest.fixture
def aligner():
    return AlignUtils()


@pytest.fixture
def input_in_fasta(path_to_test_data):
    """Returns the input data, which contains set of sequences in FASTA format for PEPTIDE1 and CHEM1"""

    with open(path_to_test_data + "input.json") as file:
        input_dict = json.load(file)

    return input_dict


@pytest.fixture
def path_to_matrices(path_to_test_data):
    """Returns paths to ROCS and monomers_map files"""

    path_to_subst_matrix = path_to_test_data + "ROCS"
    path_to_monomer_table = path_to_test_data + "monomers_map.txt"

    return [path_to_subst_matrix, path_to_monomer_table]


"""
Test the correct initiation of the AlignUtils class
"""


def test_init(aligner):

    assert aligner._AlignUtils__chars == []
    assert aligner._AlignUtils__d_enc == {}
    assert aligner._AlignUtils__d_dec == {}
    assert aligner._AlignUtils__natAA == ["G", "A", "V", "L", "I", "M", "P", "F", "W", "S", "T", "N", "Q", "Y", "C", "K", "R", "H", "D", "E"]


"""
Below are tests for __repl() function
"""


def test___repl(aligner):
    """Test replacement of non-natural AAs by unicode characters"""

    sequence = "[ClAc]FSV[Sar][Ahp]RR[NMeF]VA[Ahp]S[Bip]CGGK[NH2]"
    nnAA = re.finditer(r'\[[^\]]*\]', sequence)

    for acid in nnAA:
        aligner._AlignUtils__repl(acid.group())

    expected_chars = ['\x01', '\x02', '\x03', '\x04', '\x05', '\x06', '\x07', '\x08', '\t', '\x0b', '\x0c', '\x0e', '\x0f', '\x10', '\x11', '\x12', '\x13', '\x14', '\x15', '\x16', '\x17', '\x18', '\x19', '\x1a', '\x1b', '\x1c', '\x1d', '\x1e', '\x1f', '!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '?', '@', 'B', 'J', 'O', 'U', 'X', 'Z', '[', '\\', ']', '^', '_', '`', 'b', 'j', 'o', 'u', 'x', 'z', '{', '|', '}', '~', '\x7f', '\x80', '\x81', '\x82', '\x83', '\x84', '\x85', '\x86', '\x87', '\x88', '\x89', '\x8a', '\x8b', '\x8c', '\x8d', '\x8e', '\x8f', '\x90', '\x91', '\x92', '\x93', '\x94', '\x95', '\x96', '\x97', '\x98', '\x99', '\x9a', '\x9b', '\x9c', '\x9d', '\x9e', '\x9f', '\xa0', '¡', '¢', '£', '¤', '¥', '¦', '§', '¨', '©', 'ª', '«', '¬', '\xad', '®', '¯', '°', '±', '²', '³', '´', 'µ', '¶', '·', '¸', '¹', 'º', '»', '¼', '½', '¾', '¿', 'À', 'Á', 'Â', 'Ã', 'Ä', 'Å', 'Æ', 'Ç', 'È', 'É', 'Ê', 'Ë', 'Ì', 'Í', 'Î', 'Ï', 'Ð', 'Ñ', 'Ò', 'Ó', 'Ô', 'Õ', 'Ö', '×', 'Ø', 'Ù', 'Ú', 'Û', 'Ü', 'Ý', 'Þ', 'ß', 'à', 'á', 'â', 'ã', 'ä', 'å', 'æ', 'ç', 'è', 'é', 'ê', 'ë', 'ì', 'í', 'î', 'ï', 'ð', 'ñ', 'ò', 'ó', 'ô', 'õ', 'ö', '÷', 'ø', 'ù', 'ú', 'û', 'ü', 'ý', 'þ', 'ÿ']
    expected_d_enc = {'[ClAc]': '\x01', '[Sar]': '\x02', '[Ahp]': '\x03', '[NMeF]': '\x04', '[Bip]': '\x05', '[NH2]': '\x06'}
    expected_d_dec = {'\x01': '[ClAc]', '\x02': '[Sar]', '\x03': '[Ahp]', '\x04': '[NMeF]', '\x05': '[Bip]', '\x06': '[NH2]'}

    assert aligner._AlignUtils__chars == expected_chars
    assert aligner._AlignUtils__d_enc == expected_d_enc
    assert aligner._AlignUtils__d_dec == expected_d_dec


def test___repl_no_more_letters(aligner, path_to_test_data):
    """Test the raise of exception, when there is no more place for additional symbols"""

    # The file contains more monomers of nnAAs (most of which are generated artificially) that the capacity of the encoder
    path_to_artificial_seqs = path_to_test_data + "artificial_seqs.json"

    with open(path_to_artificial_seqs) as file:

        # Read the file
        artificial_seqs = json.load(file)
        artificial_seqs = artificial_seqs["artificial_seqs"]

        with pytest.raises(Exception) as excinfo:

            # Iterate through the sequences
            for entry in artificial_seqs:

                # Remove dots from HELM seqs
                entry = entry.replace(".", "")

                # Find all the nnAA monomers in the sequence
                findings = re.finditer(r'\[[^\]]*\]', entry)

                # Encode the nnAA monomers
                for acid in findings:
                    aligner._AlignUtils__repl(acid)

        assert str(excinfo.value) == "There are no more letters in mapping array"


"""
Below are tests for encode_alignment_sequences() function
"""


def test_encode_alignment_sequences(aligner, path_to_test_data, input_in_fasta):
    """ Test the correct encoding of the nnAAs """

    for key in input_in_fasta:
        data = input_in_fasta[key]
        encoded_file = aligner.encode_alignment_sequences(data, key, 'TEST')
        test_file = path_to_test_data + f"{key}_TEST_mafft.txt"
        with open(encoded_file) as output_file, open(test_file) as expected_output:
            output_file = output_file.readlines()
            expected_output = expected_output.readlines()
            i = 0
            for output_line, expected_line in zip(output_file, expected_output):
                i += 1
                assert output_line == expected_line

        aligner.clear_symbols()
        os.remove(encoded_file)


@pytest.mark.xfail
def test_encode_alignment_sequences_fails_no_cleaning(aligner, path_to_test_data):
    """ Test the correct encoding of the nnAAs """

    with open(path_to_test_data + "input.json") as file:
        input_dict = json.load(file)

    for key in input_dict:
        data = input_dict[key]
        encoded_file = aligner.encode_alignment_sequences(data, key, 'TEST')
        test_file = path_to_test_data + f"{key}_TEST_mafft.txt"
        with open(encoded_file) as output_file, open(test_file) as expected_output:
            output_file = output_file.readlines()
            expected_output = expected_output.readlines()

            # There will be AssertionError in this test, so we need to remove the encoded file beforehand
            os.remove(encoded_file)

            i = 0
            for output_line, expected_line in zip(output_file, expected_output):
                i += 1
                assert output_line == expected_line

        # Symbols are not cleared, so aligner will use another set of symbols, to encode nnAAs, which will result in AssertionError


def test_encode_alignment_sequences_fails_encoder_overfilled(aligner, path_to_test_data):
    """Test that error is raised due to the lack of letters for encoding"""

    path_to_artificial_seqs = path_to_test_data + "artificial_seqs.json"

    with open(path_to_artificial_seqs) as file:
        artificial_seqs = json.load(file)
        artificial_seqs = artificial_seqs["artificial_seqs"]

        with pytest.raises(Exception) as excinfo:
            aligner.encode_alignment_sequences(artificial_seqs, "PEPTIDE", "TEST")

        assert str(excinfo.value) == "There are no more letters in mapping array"


"""
Below are tests for decode_mafft() function
"""


def test_decode_mafft(aligner, path_to_test_data):
    """
    encoded = ['> PEPTIDE1\n\x01--------------\x02RRRRCPLYIS\x03DPVCRRRR\x04\n', '> PEPTIDE1\n\x01-------------------\x05PLYISYDPV\x06----\x04\n', '> PEPTIDE1\n\x01--------------\x02RRRR\x07PLYISYDPV\x08RRRR\x04\n', '> PEPTIDE1\n\x06----------------RRRCPLYISYDPVCRRR\x06\x04\n', '> PEPTIDE1\n\x01RQIKIWFQNRRMKWKKG\x02\t\x0bPLYISYDPVC---R\x04\n', '> PEPTIDE1\n\x01--------------\x02RRRR\x0b\x0cLYISYDPVCRRRR\x04\n', '> PEPTIDE1\n\x01----------------\x02RR\x0bPLYISYDPV\x0e--RR\x04\n', '> PEPTIDE1\n\x0f------FSV-------\x10\x11RR\x12VA\x11S\x13CGG----K\x04\n', '> PEPTIDE1\n--------------------ACAKCA\t---------\n', '> PEPTIDE1\n\x01--------------K\x14\x14\x14\x14\x0b\x15LYI\x16YDPVC\x14\x14\x14\x14\x04\n', '> PEPTIDE1\nH--\x17QGTFTSEYSKYLDERAAKDFVQWLLN----K\x04\n', '> PEPTIDE1\n\x01-------------------\x0b\x18LFI\x05YDPVC---\x19\x04\n', '> PEPTIDE1\n\x1a------F------------R\x1bLY\x1c\x12----------\n', '> PEPTIDE1\n\x1d-------------------PKLY\x1e\x12--------\x1f\x04\n', '> PEPTIDE1\n\x1d-------------------PKLY\x1e\x12--------!\x04\n']
    expected_output = ['> PEPTIDE1\n[Ac]--------------[LysN3]RRRRCPLYIS[NMeY]DPVCRRRR[NH2]\n', '> PEPTIDE1\n[Ac]-------------------[dApe]PLYISYDPV[Ape]----[NH2]\n', '> PEPTIDE1\n[Ac]--------------[LysN3]RRRR[dE]PLYISYDPV[Dab]RRRR[NH2]\n', '> PEPTIDE1\n[Ape]----------------RRRCPLYISYDPVCRRR[Ape][NH2]\n', '> PEPTIDE1\n[Ac]RQIKIWFQNRRMKWKKG[LysN3][AEEA][dC]PLYISYDPVC---R[NH2]\n', '> PEPTIDE1\n[Ac]--------------[LysN3]RRRR[dC][Prot3OH]LYISYDPVCRRRR[NH2]\n', '> PEPTIDE1\n[Ac]----------------[LysN3]RR[dC]PLYISYDPV[Pen]--RR[NH2]\n', '> PEPTIDE1\n[ClAc]------FSV-------[Sar][Ahp]RR[NMeF]VA[Ahp]S[Bip]CGG----K[NH2]\n', '> PEPTIDE1\n--------------------ACAKCA[AEEA]---------\n', '> PEPTIDE1\n[Ac]--------------K[dR][dR][dR][dR][dC][dProt3Ph]LYI[aMeS]YDPVC[dR][dR][dR][dR][NH2]\n', '> PEPTIDE1\nH--[dS]QGTFTSEYSKYLDERAAKDFVQWLLN----K[NH2]\n', '> PEPTIDE1\n[Ac]-------------------[dC][Prot3Ph]LFI[dApe]YDPVC---[Hag][NH2]\n', '> PEPTIDE1\n[2moPyr]------F------------R[4Pal]LY[Nva][NMeF]----------\n', '> PEPTIDE1\n[dF]-------------------PKLY[Nle][NMeF]--------[hE][NH2]\n', '> PEPTIDE1\n[dF]-------------------PKLY[Nle][NMeF]--------[Apm][NH2]\n']

    # First, we need to create encoding/decoding dictionary. We'll use encode_alignment_sequences() function for this
    encoded_file = aligner.encode_alignment_sequences(expected_output, 'PEPTIDE1', 'TEST')

    # encode_alignment_sequences() creates file, so we need to rmove it
    os.remove(encoded_file)

    # Decode aligned sequences
    output = aligner.decode_mafft(encoded)

    assert output == expected_output
    """
    path = path_to_test_data + "decode.json"

    with open(path, encoding="utf-8") as file:
        data = json.load(file)

        for key in data:
            encoded = data[key]["input"]
            expected_output = data[key]["expected_output"]

            # First, we need to create encoding/decoding dictionary. We'll use encode_alignment_sequences() function for this
            encoded_file = aligner.encode_alignment_sequences(expected_output, key, 'TEST')

            # encode_alignment_sequences() creates file, so we need to remove it
            os.remove(encoded_file)

            # Decode aligned sequences
            output = aligner.decode_mafft(encoded)
            aligner.clear_symbols()

            assert output == expected_output


@pytest.mark.xfail
def test_decode_mafft_fails_no_cleaning(aligner, path_to_test_data):
    """Test that aligner will use different decoding symbols, if they are not cleared"""
    path = path_to_test_data + "decode.json"

    with open(path, encoding="utf-8") as file:
        data = json.load(file)

        for key in data:
            encoded = data[key]["input"]
            expected_output = data[key]["expected_output"]

            # First, we need to create encoding/decoding dictionary. We'll use encode_alignment_sequences() function for this
            encoded_file = aligner.encode_alignment_sequences(expected_output, key, 'TEST')

            # encode_alignment_sequences() creates file, so we need to remove it
            os.remove(encoded_file)

            # Decode aligned sequences
            output = aligner.decode_mafft(encoded)

            # Here we do not perform clearing of the symbols, which will result in assertion error

            assert output == expected_output


"""
Below are tests for clear_symbols() function
"""


def test_clear_symbols(aligner):

    # assign values to the private variables
    aligner._AlignUtils__chars = ['\x01', '\x02', '\x03']
    aligner._AlignUtils__d_enc = {'[Ac]': '\x01', '[LysN3]': '\x02', '[NMeY]': '\x03'}
    aligner._AlignUtils__d_dec = {'\x01': '[Ac]', '\x02': '[LysN3]', '\x03': '[NMeY]'}

    # Clear the contents of the variables
    aligner.clear_symbols()

    assert aligner._AlignUtils__chars == []
    assert aligner._AlignUtils__d_enc == {}
    assert aligner._AlignUtils__d_dec == {}


"""
Below are tests for get_unicode_char() function
"""


def test_get_unicode_char():
    """Test the correct conversion of strings to unicode characters"""

    inp = ["0100", "020E", "084C", "0391", "0405", "050F", "062E", "071F"]
    expected_output = ["Ā", "Ȏ", "ࡌ", "Α", "Ѕ", "ԏ", "خ", "ܟ"]

    for input, expected in zip(inp, expected_output):
        output = AlignUtils.get_unicode_char(input)
        assert output == expected


"""
Below are tests for read_fasta() function
"""


def test_read_fasta(aligner):
    """ Test reading of the FASTA input"""

    input_array = ["> PEPTIDE1\n[ClAc]FSV[Sar][Ahp]RR[NMeF]VA[Ahp]S[Bip]CGGK[NH2]\n",
                   "> PEPTIDE1\nACAKCA\n",
                   "> PEPTIDE1\n[TCO][PEG6NH2acid]\n",
                   "> PEPTIDE1\n[C18diacid][gE][AEEA][AEEA]\n"]

    for ind, output in enumerate(aligner.read_fasta(input_array)):
        expected_output = tuple(input_array[ind].strip().split("\n"))
        assert output == expected_output


"""
Below are tests for create_substitution_matrix() function
"""


def test_create_substitution_matrix(aligner, input_in_fasta, path_to_test_data, path_to_matrices):
    """Test the creation of substitution matrix"""

    # Assign paths to the ROCS and monomers_map files
    path_to_subst_matrix, path_to_monomer_table = path_to_matrices

    # Iterate through input data (for PEPTIDE1 and CHEM1)
    for key in input_in_fasta:

        # Read the data
        data = input_in_fasta[key]

        # Artificial timestamp
        stamp = "TEST"

        # Encode the data
        encoded_file = aligner.encode_alignment_sequences(data, key, stamp)

        # Create substitution matrix
        subst_matrix_file = aligner.create_substitution_matrix(path_to_subst_matrix,
                                                               path_to_monomer_table,
                                                               key,
                                                               stamp)
        # Name of the file for comparison of matrices
        expected_file = f"{key}_EXPECTED_matrix.txt"
        with open(subst_matrix_file) as output, open(path_to_test_data + expected_file) as expected_file:
            output_data = output.readlines()
            expected_data = expected_file.readlines()

            for output_line, expected_line in zip(output_data, expected_data):
                assert output_line == expected_line

        # Clear symbols of the aligner
        aligner.clear_symbols()

        # Remove file with encoded data and substitution matrix file
        os.remove(encoded_file)
        os.remove(subst_matrix_file)


"""
Below are tests for run_mafft_utility() function
"""


def test_run_mafft_utility(aligner, path_to_test_data, input_in_fasta, path_to_matrices):
    """Test running of the MAFFT program"""

    # Assign paths to the ROCS and monomers_map files
    path_to_subst_matrix, path_to_monomer_table = path_to_matrices

    # Change  to "alignment/mafft/bin/ginsi" for develop
    mafft_binary = "ginsi"
    key = "PEPTIDE1"
    stamp = "TEST"

    data = input_in_fasta[key]
    encoded_data = aligner.encode_alignment_sequences(data, key, stamp)
    subst_matrix = aligner.create_substitution_matrix(path_to_subst_matrix, path_to_monomer_table, key, stamp)

    mafft_output, mafft_error = aligner.run_mafft_utility(encoded_data, mafft_binary=mafft_binary, matrix_file=subst_matrix, realign=False,
                                                          gap_opening_penalty=1, gap_extension_penalty=0)

    mafft_output = mafft_output.split("\n")

    with open(path_to_test_data + "PEPTIDE1_mafft_aligned.txt") as expected_file:
        expected_data = expected_file.readlines()

        for output_line, expected_line in zip(mafft_output, expected_data):
            expected_line = expected_line.rstrip("\n")
            assert output_line == expected_line

    os.remove(encoded_data)
    os.remove(subst_matrix)


"""
Below are tests for parse_naa_in_fasta() function
"""


def test_parse_naa_in_fasta(aligner):
    """Test the correct parsing of sequence with non-natural AAs"""

    seq = "[ClAc]FSV[Sar][Ahp]RR[NMeF]VA[Ahp]S[Bip]CGGK[NH2]"
    expected_list = ["[ClAc]", "[Sar]", "[Ahp]", "[NMeF]", "[Ahp]", "[Bip]", "[NH2]"]

    out_list = aligner.parse_naa_in_fasta(seq)

    assert expected_list == out_list


def test_parse_naa_in_fasta_smiles(aligner):
    """Test the correct parsing of sequence with nested non-natural AAs"""

    seq = "[FLac][NMeA]KK[cBuA]GY[[R2]C(=O)[C@H](CCC1CC1)N[R1]][NMeF][amBCP]"
    expected_list = ["[FLac]", "[NMeA]", "[cBuA]", "[[R2]C(=O)[C@H](CCC1CC1)N[R1]]", "[NMeF]", "[amBCP]"]

    out_list = aligner.parse_naa_in_fasta(seq)

    assert expected_list == out_list


def test_parse_naa_in_fasta_empty(aligner):
    """Test the correct parsing of sequence without non-natural AAs"""

    seq = "ABCDEFGHIJ"
    expected_list = []

    out_list = aligner.parse_naa_in_fasta(seq)

    assert expected_list == out_list


def test_parse_naa_in_fasta_fails(aligner):
    """Test the raising of the error, in case of invalid notation (inconsistent number of opening and closing brackets)"""

    seq = "[FLac][NMeAKK[cBuA]GY[[R2]C(=O)[C@H](CCC1CC1)N[R1]][NMeF][amBCP]"

    with pytest.raises(Exception) as excinfo:
        aligner.parse_naa_in_fasta(seq)

    assert str(excinfo.value) == f"Non-valid notation of non-natural amino acids in the following sequence {seq}"
