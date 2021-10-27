import os
import json
import glob
import pytest
from datetime import datetime

from alignment.AlignSubPeptides import align_sub_peptides, split_sub_peptides, extract_sub_peptide, \
    helm2fasta, fasta2helm, convert2helm, get_common_subpeptides, get_aligned_sequences, get_alignment_score, \
    cleanup


# Pytest fixtures
@pytest.fixture
def path_to_data():
    # Return path to data folder

    path = "tests/unit_tests/resources/test_AlignSubPeptides/"
    return path


@pytest.fixture
def helm_sequences():
    """Returns prefixes of files with raw HELM sequences"""

    return ['example_', 'for_alignment1_', 'for_alignment2_', 'for_alignment3_', 'for_alignment_without_nh2_',
            'helmtest_', 'input_data_', 'multi_chain_peps_', 'multi_chem_pep_', "test_data_"]


@pytest.fixture
def non_empty_dict():
    """Returns non-empty HELM dictionary"""

    return {"PEPTIDE1": "PEPTIDE1{[c3amCb1c].[d1Nal].R.K.[Nle].Y.[Nle].[NMeF]}$$$$"}


@pytest.fixture
def data_for_alignment(path_to_data):
    """
    Return data, needed to test the correct alignment
      - input data, expected output (alignment and score), matrix and monomers table"""
    input_path = path_to_data + "align_input.txt"
    output_path = path_to_data + "align_output.json"

    matrix = path_to_data + "ROCS"
    monomers = path_to_data + "monomers_map.txt"

    with open(input_path)as inp_file, open(output_path) as out_file:
        input_data = inp_file.readlines()
        data = ""
        for line in input_data:
            data += line

        # Expected output is generated with the help of ginsi method
        expected_alignment = json.load(out_file)
        expected_score = {"PEPTIDE1": -1.1537037037037037, "CHEM1": -19.375, "PEPTIDE2": -7.532608695652174,
                          "CHEM2": 0.0, "CHEM3": None}

        return data, expected_alignment, expected_score, matrix, monomers


@pytest.fixture
def data_for_realignment(path_to_data):
    """Return data, needed to test the correct realignment
    - input data, expected output (alignment and score), matrix and monomers table"""

    aligned_path = path_to_data + "align_output.json"
    new_seqs = path_to_data + "new_sequences.txt"
    realign_path = path_to_data + "realign_output.json"

    matrix = path_to_data + "ROCS"
    monomers = path_to_data + "monomers_map.txt"

    with open(aligned_path)as align, open(realign_path) as realign, open(new_seqs) as new_file:
        aligned_seqs = json.load(align)
        aligned_data = ""
        for key in aligned_seqs:
            aligned_data += aligned_seqs[key]

        new_sequences = new_file.readlines()
        new_data = ""
        for line in new_sequences:
            new_data += line

        # Expected output is generated with the help of ginsi method
        expected_alignment = json.load(realign)
        expected_score = {"PEPTIDE1": 216.61280193236715,
                          "PEPTIDE2": 276.5110946745562}

        return aligned_data.strip("\n"), new_data, expected_alignment, expected_score, matrix, monomers


"""
Below are test for align_sub_peptides() function
"""


def test_align_sub_peptides(path_to_data, data_for_alignment):
    """Test the correct alignment"""

    data, expected_alignment, expected_score, matrix, monomers = data_for_alignment

    output, err, score = align_sub_peptides(data, gap_opening_penalty=1, gap_extension_penalty=0, polymer_to_align=None, path_to_mafft="ginsi",
                                            path_to_subst_matrix=matrix, path_to_monomer_table=monomers)

    # Remove temporary files
    files_to_remove = []
    for key in output:
        files = glob.glob(f'*{key}*')
        files_to_remove.extend(files)

    for file in files_to_remove:
        os.remove(file)

    # Assert the correctness of the output
    for key in output:
        assert output[key] == expected_alignment[key]
        assert score[key] == expected_score[key]


def test_align_sub_peptides_single_polymer(data_for_alignment):
    """Test the correct alignment"""

    data, expected_alignment, expected_score, matrix, monomers = data_for_alignment

    # Call align_sub_peptides() and specify polymer_to_align, which is not present in the input data
    output, err, score = align_sub_peptides(data, gap_opening_penalty=1, gap_extension_penalty=0, polymer_to_align="PEPTIDE1", path_to_mafft="ginsi",
                                            path_to_subst_matrix=matrix, path_to_monomer_table=monomers)

    # Remove temporary files
    files_to_remove = glob.glob(f'*{"PEPTIDE1"}*')

    for file in files_to_remove:
        os.remove(file)

    assert len(output) == 1
    assert len(score) == 1
    assert output["PEPTIDE1"] == expected_alignment["PEPTIDE1"]
    assert score["PEPTIDE1"] == expected_score["PEPTIDE1"]


def test_align_sub_peptides_missing_polymer(data_for_alignment):
    """Test the correct alignment"""

    data, expected_alignment, expected_score, matrix, monomers = data_for_alignment

    # Call align_sub_peptides() and specify polymer_to_align, which is not present in the input data
    output, err, score = align_sub_peptides(data, gap_opening_penalty=1, gap_extension_penalty=0, polymer_to_align="PEPTIDE3", path_to_mafft="ginsi",
                                            path_to_subst_matrix=matrix, path_to_monomer_table=monomers)

    assert output is None
    assert err is None
    assert score is None


def test_align_sub_peptides_realignment(path_to_data, data_for_realignment):
    """Test the correct realignment"""

    aligned_data, new_data, expected_alignment, expected_score, matrix, monomers = data_for_realignment
    output, err, score = align_sub_peptides(aligned_data, new_data, gap_opening_penalty=1, gap_extension_penalty=0,
                                            polymer_to_align=None, path_to_mafft="ginsi",
                                            path_to_subst_matrix=matrix, realign_method="add",
                                            path_to_monomer_table=monomers)

    # Remove temporary files
    files_to_remove = []
    for key in output:
        files = glob.glob(f'*{key}*')
        files_to_remove.extend(files)

    for file in files_to_remove:
        os.remove(file)

    # Assert, that for a given output, 2 common polymers were aligned
    assert len(output) == 2

    # Assert the correctness of the output
    for key in output:
        assert output[key] == expected_alignment[key]
        assert score[key] == expected_score[key]


def test_align_sub_peptides_realignment_specified_polymer_exists(path_to_data, data_for_realignment):
    """Test the correct realignment, with specified polymer"""

    polymer_to_align = "PEPTIDE2"
    aligned_data, new_data, expected_alignment, expected_score, matrix, monomers = data_for_realignment
    output, err, score = align_sub_peptides(aligned_data, new_data, gap_opening_penalty=1, gap_extension_penalty=0,
                                            polymer_to_align=polymer_to_align, path_to_mafft="ginsi",
                                            path_to_subst_matrix=matrix, realign_method="add",
                                            path_to_monomer_table=monomers)

    # Remove temporary files
    files_to_remove = []
    for key in output:
        files = glob.glob(f'*{key}*')
        files_to_remove.extend(files)

    for file in files_to_remove:
        os.remove(file)

    # Assert, that for a given output, 2 common polymers were aligned
    assert len(output) == 1

    # Assert the correctness of the output
    assert output[polymer_to_align] == expected_alignment[polymer_to_align]
    assert score[polymer_to_align] == expected_score[polymer_to_align]


def test_align_sub_peptides_realignment_specified_polymer_missing(path_to_data, data_for_realignment):
    """Test that no alignment is performed, if the specified polymer is missing"""

    # PEPTIDE3 is not present in either input data
    polymer_to_align = "PEPTIDE3"
    aligned_data, new_data, expected_alignment, expected_score, matrix, monomers = data_for_realignment
    output, err, score = align_sub_peptides(aligned_data, new_data, gap_opening_penalty=1, gap_extension_penalty=0,
                                            polymer_to_align=polymer_to_align, path_to_mafft="ginsi",
                                            path_to_subst_matrix=matrix, realign_method="add",
                                            path_to_monomer_table=monomers)

    assert output is None
    assert err is None
    assert score is None


"""
Below are tests for function split_sub_peptides()
"""


def test_split_sub_peptides(helm_sequences, path_to_data):

    """This function goes through different test input files.
       Path to folder with data is stored in dir variable, files prefixes are stored in helm_sequences() fixture
            - Input is stored in 'prefix_input.txt'
            - Output is stored in 'prefix_*.json', where * regexp specifies the type os subpeptide used for splitting
    """

    dir = path_to_data + "test_split_sub_peptides"

    # List of files in directory
    files = os.listdir(dir)

    # This loop iterates over prefixes of input files
    for prefix in helm_sequences:

        # List of files output files
        output_files = []

        # List of respective subpeptides
        subpeptides = []

        for file in files:
            if prefix + "input" in file:

                # Input file
                input_file = file
            elif prefix + "all" in file:

                # Output file with all subpeptides
                output_with_all_subpeptides = file
            elif prefix in file:

                # Extract peptides
                peptide = file.split(prefix)[1].split("_")[0]
                subpeptides.append(peptide)

                output_files.append(file)

        # Add None, to the list of peptides, to check the output with None parameter on output file, which contains all the peptides
        subpeptides = [None, *subpeptides]
        path_in = dir + "/" + input_file
        with open(path_in) as input_lines:
            input = ""
            for line in input_lines:
                input += line

            ind = 0
            for peptide in subpeptides:

                # If None, then pass to the function input data and None, and check on output file, which contains all the peptides
                if peptide is None:
                    result = split_sub_peptides(input, None)
                    file = output_with_all_subpeptides
                else:
                    result = split_sub_peptides(input, peptide)
                    file = output_files[ind]
                    ind += 1
                path_out = dir + "/" + file

                with open(path_out) as output:
                    expected_output = json.load(output)
                    assert result == expected_output


@pytest.mark.xfail
def test_split_sub_peptides_fails():
    """Test with wrong input"""
    split_sub_peptides("Some Text")


"""
Below are tests for function extract_sub_peptide()
"""


def test_extract_subpeptide_empty_dict():
    """Adds sequence in empty dictionary"""
    helm_dict = {}
    helm_dict = extract_sub_peptide("PEPTIDE1", helm_dict, "[ClAc].F.R.Y.L.Y.[Ahp].F.C.G.K.K.[NH2]")

    # Supposed output
    result = {"PEPTIDE1": "PEPTIDE1{[ClAc].F.R.Y.L.Y.[Ahp].F.C.G.K.K.[NH2]}$$$$"}

    assert helm_dict == result


def test_extract_sub_peptide_nonempty_dict(non_empty_dict):
    """Adds sequence to the existing subpeptide in the dictionary"""

    helm_dict = extract_sub_peptide("PEPTIDE1", non_empty_dict, "[ClAc].F.R.Y.L.Y.[Ahp].F.C.G.K.K.[NH2]")

    # Supposed output
    result = {"PEPTIDE1": "PEPTIDE1{[c3amCb1c].[d1Nal].R.K.[Nle].Y.[Nle].[NMeF]}$$$$\nPEPTIDE1{[ClAc].F.R.Y.L.Y.[Ahp].F.C.G.K.K.[NH2]}$$$$"}

    assert helm_dict == result


def test_extract_sub_peptide_add_new_subpeptide(non_empty_dict):
    """Adds new subpeptide to the HELM dictionary"""

    helm_dict = extract_sub_peptide("PEPTIDE2", non_empty_dict, "[C18diacid].[gE].[AEEA].[AEEA]")

    # Supposed output
    result = {"PEPTIDE1": "PEPTIDE1{[c3amCb1c].[d1Nal].R.K.[Nle].Y.[Nle].[NMeF]}$$$$",
              "PEPTIDE2": "PEPTIDE2{[C18diacid].[gE].[AEEA].[AEEA]}$$$$"}

    assert helm_dict == result


@pytest.mark.xfail
def test_extract_sub_peptide_fails():
    """Test with wrong input"""

    extract_sub_peptide("Some Text")


"""
Below are tests for function helm2fasta()
"""


def test_helm2fasta():
    """Converts sequence from HELM format into FASTA"""

    # Input HELM string
    helm = "PEPTIDE1{H.[dS].Q.G.T.F.T.S.E.Y.S.K.Y.L.D.E.R.A.A.K.D.F.V.Q.W.L.L.N.K.[NH2]}|PEPTIDE2{[C18diacid].[gE].[AEEA].[AEEA]}$PEPTIDE1,PEPTIDE2,12:R3-4:R2|PEPTIDE1,PEPTIDE1,16:R3-20:R3$$$V2.0"

    # Supposed output
    result = "> PEPTIDE1\nH[dS]QGTFTSEYSKYLDERAAKDFVQWLLNK[NH2]\n> PEPTIDE2\n[C18diacid][gE][AEEA][AEEA]\n"

    fasta = helm2fasta(helm)

    assert fasta == result


def test_helm2fasta_wrong_helm():
    """Test with wrong input"""

    with pytest.raises(Exception) as excinfo:
        helm2fasta('Just Text')

    assert str(excinfo.value) == "Invalid HELM notation"


"""
Below are tests for function convert2helm()
"""


def test_convert2helm():
    """Tests transformation from FASTA format into HELM"""

    # Input FASTA data
    fasta = '> PEPTIDE1\nH[dS]QGTFTSEYSKYLDERAAKDFVQWLLNK[NH2]\n> PEPTIDE1\n[2moPyr]FR[4Pal]LY[Nva][NMeF]'

    # Supposed output
    result = ["PEPTIDE1{H.[dS].Q.G.T.F.T.S.E.Y.S.K.Y.L.D.E.R.A.A.K.D.F.V.Q.W.L.L.N.K.[NH2]}$$$$", "PEPTIDE1{[2moPyr].F.R.[4Pal].L.Y.[Nva].[NMeF]}$$$$"]
    helm = convert2helm(fasta)

    assert helm == result


# Should we leave this?
@pytest.mark.xfail
def test_convert2helm_failed():
    """Test with wrong input"""

    # The following line doesn't lead to error. Maybe input validation is needed.
    convert2helm("Some Text")


"""
Below are tests for function fasta2helm()
"""


def test_fasta2helm():
    """Tests transformation of simple FASTA string into HELM string"""

    # HELM sequence in FASTA format
    helm = "H[dS]QGTFTSEYSKYLDERAAKDFVQWLLNK[NH2]"

    # HELM sequence which is supposed to be returned by fasta2helm() function
    output = "H.[dS].Q.G.T.F.T.S.E.Y.S.K.Y.L.D.E.R.A.A.K.D.F.V.Q.W.L.L.N.K.[NH2]"
    fasta = fasta2helm(helm)

    assert fasta == output


def test_fasta2helm_corrupted_sequence():
    """Tests that corrupted sequence results in error"""

    # corrupted HELM sequence
    helm = "H[dSQGTFTSEYSKYLDERAAKDFVQWLLNK"

    with pytest.raises(Exception) as excinfo:
        fasta2helm(helm)

    expected_message = f"Non-valid notation of non-natural amino acids in the following sequence {helm}"

    assert str(excinfo.value) == expected_message


"""
Below are tests for function get_common_subpeptides()
"""


def test_get_common_subpeptides():
    """Test that get_common_subpeptides() function returns intersection of 2 lists"""

    # Lists of peptides
    l1 = ['PEPTIDE1', 'PEPTIDE2']
    l2 = ['PEPTIDE2', 'PEPTIDE3']
    com_list = get_common_subpeptides(l1, l2)

    assert com_list == ['PEPTIDE2']


"""
Below are tests for function get_aligned_sequences()
"""


def test_get_aligned_sequences():
    """Tests the splitting of the whole MAFFT output into the separate entities"""

    # Output of the MAFFT program
    aligned = "> PEPTIDE1\n[Ac]--------------[LysN3]RRRR[dE]PLYISYDPV[Dab]RRRR[NH2]\n> PEPTIDE1\n[Ac]--------------[LysN3]RRRR[dC][Prot3OH]LYISYDPVCRRRR[NH2]\n"

    # Supposed output of the get_aligned_sequences()
    result = ["> PEPTIDE1\n[Ac]--------------[LysN3]RRRR[dE]PLYISYDPV[Dab]RRRR[NH2]\n",
              "> PEPTIDE1\n[Ac]--------------[LysN3]RRRR[dC][Prot3OH]LYISYDPVCRRRR[NH2]\n"]

    output = get_aligned_sequences(aligned)

    assert result == output


"""
Below are tests for function cleanup()
"""


def test_cleanup():
    """Test cleanup() function"""

    # Create timestamp for file
    timestamp = datetime.now().isoformat(timespec='milliseconds')

    # Create files with timestamp
    with open("test_cleanup_file_1_" + timestamp + "_.txt", 'w') as file1,\
            open("test_cleanup_file_2_" + timestamp + "_.txt", 'w') as file2:
        file1.write('Some text')
        file2.write('Some other text')
        pass

    # Delete these files
    cleanup(timestamp)

    deleted = True
    files = os.listdir()

    # Iterate through the files in directory, and assign False to "deleted" variable, if file with the timestamp exists
    for name in files:
        if timestamp in name:
            deleted = False
            break

    assert deleted is True


"""
Below are tests for function get_alignment_score()
"""


def test_get_alignment_score(path_to_data):
    """Test scoring of the alignment"""

    path_to_test_data = path_to_data + "for_scoring.txt"
    path_to_matrix = path_to_data + "PEPTIDE1_SCORE_matrix.txt"

    with open(path_to_test_data) as data:
        data = "".join(data.readlines())
        expected_score = -1.1537037037037037

        result = get_alignment_score(data, path_to_matrix)

        assert result == expected_score


def test_get_alignment_score_empty(path_to_data):
    """Test that on attempt to score alignment, consisting of one sequence, the function returns None"""

    path_to_test_data = path_to_data + "for_scoring.txt"
    path_to_matrix = path_to_data + "PEPTIDE1_SCORE_matrix.txt"

    with open(path_to_test_data) as data:
        data = "".join(data.readlines()[:2])

        result = get_alignment_score(data, path_to_matrix)

        assert result is None
