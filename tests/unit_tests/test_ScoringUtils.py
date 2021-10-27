import pytest

from pymsa import MSA, SumOfPairs
from alignment.scoring.ScoringUtils import CustomMatrix, Scoring


@pytest.fixture
def path_to_resources():
    # Return path to the test data

    path = "tests/unit_tests/resources/test_ScoringUtils/"
    return path


@pytest.fixture
def path_to_matrix(path_to_resources):
    # Return path to the substitution matrix

    path = path_to_resources + "PEPTIDE1_FOR_SCORING_matrix.txt"
    return path


@pytest.fixture
def custom_matrix(path_to_matrix):
    # Return custom matrix

    matrix = CustomMatrix(path_to_matrix)
    return matrix


@pytest.fixture
def expected_matrix(path_to_resources):
    # Return expected distance matrix

    path = path_to_resources + "expected_matrix.txt"
    with open(path, encoding="utf-8") as file:
        file = file.readlines()
        matrix = dict()
        for line in file:
            k1, k2, value = line.split("\t")
            matrix[(chr(int(k1)), chr(int(k2)))] = float(value)
        return matrix


@pytest.fixture
def mafft_output(path_to_resources):
    # return aligned sequences

    path = path_to_resources + "aligned_seqs.txt"

    with open(path) as file:
        aligned_seqs = "".join(file.readlines())
        return aligned_seqs


"""
Below are tests for CustomMatrix class
"""


def test_CustomMatrix_init(custom_matrix, expected_matrix):

    assert len(custom_matrix.distance_matrix) == len(expected_matrix)

    for key in custom_matrix.distance_matrix:
        assert custom_matrix.distance_matrix[key] == expected_matrix[key]


"""
Below are tests for read_matrix_from_file() method
"""


def test_read_matrix_from_file(path_to_matrix, expected_matrix):

    matrix = CustomMatrix.read_matrix_from_file(path_to_matrix)
    assert matrix == expected_matrix


"""
Below are tests for Scoring class
"""


def test_mafft_output_to_msa(mafft_output):
    """Test the correct transformation of the MAFFT output into MSA object"""

    msa = Scoring.mafft_output_to_msa(mafft_output)

    assert msa.is_valid

    aligned_seqs = mafft_output.strip().split("\n")
    ids = [i for i in aligned_seqs[::2]]
    seqs = [i for i in aligned_seqs[1::2]]

    expected_msa = MSA(seqs, ids)

    for id_out, id_exp, seq_out, seq_exp in zip(msa.ids, expected_msa.ids, msa.sequences, expected_msa.sequences):

        assert id_out == id_exp
        assert seq_out == seq_exp


"""
Below are tests for sum_of_pairs() method
"""


def test_sum_of_pairs(path_to_matrix, mafft_output):
    """Test the correct scoring of alignment"""

    msa = Scoring.mafft_output_to_msa(mafft_output)
    score = Scoring.sum_of_pairs(msa, path_to_matrix)

    matrix = CustomMatrix(path_to_matrix)
    expected_score = SumOfPairs(msa, matrix).compute()

    assert score == expected_score
