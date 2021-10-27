import csv
import re
from pymsa import MSA, SumOfPairs
from pymsa import SubstitutionMatrix


class CustomMatrix(SubstitutionMatrix):
    """
    Custom class, which inherits functionality of SubstitutionMatrix class from pyMSA package
    Helps to create substitution matrix from file, which is created for MAFFT program

    It is much faster, as substitution matrix file contains only necessary monomers
    """
    def __init__(self, path, gap_penalty: int = -10, gap_character: str = "-"):
        """
        Default values for gap_penalty and gap_character are the same,
        as the corresponding values of substitution matrices in pyMSA package (PAM250, BLOSUM62)

        :param path:            path to the custom substitution matrix, which is generated for MAFFT program
        :param gap_penalty:     gap penalty value
        :param gap_character:   gap symbol
        :return:                distance matrix - dictionary, where keys are 2-element tuple of symbols, and values are scores of their alignment
        """
        super(CustomMatrix, self).__init__(gap_penalty, gap_character)
        self.distance_matrix = self.read_matrix_from_file(path)

    @staticmethod
    def read_matrix_from_file(path_to_matrix: str) -> dict:
        """
        Parses substitution matrix file, extracts symbols and their alignment score from each line

        :param path_to_matrix:  path to the custom substitution matrix
        :return:                distance matrix
        """
        distance_matrix = {}

        with open(path_to_matrix) as csv_file:
            reader = csv.reader(csv_file)
            for row in reader:
                tokens = re.split(r'\s+', row[0])

                # Decoding of characters uses the same approach, as get_unicode_char() method of AlignUtils
                symb1 = chr(int(tokens[0], 16))
                symb2 = chr(int(tokens[1], 16))

                """
                As substitution matrix is square it contains redundant information,
                so the following condition checks, that distance matrix will only contain unique pairs of symbols
                """
                if (symb1, symb2) not in distance_matrix and (symb2, symb1) not in distance_matrix:
                    distance_matrix.update({(symb1, symb2): int(tokens[2])})

        return distance_matrix


class Scoring:
    """
    Methods for data manipulation and calculation of scores
    """

    @staticmethod
    def mafft_output_to_msa(mafft_output):
        """
        Transforms the output of MAFFT program into the MSA object.
        It is neccessary as methods of pyMSA package work with these objects as input, rather than raw sequences

        :param mafft_output:    output of MAFFT program
        :return:                MSA object of aligned sequences
        """
        mafft_lines = mafft_output.splitlines()

        headers = [line for line in mafft_lines if line.startswith(">")]
        sequences = ["".join(sequence.split("\n")[1:]) for sequence in mafft_output.split(">") if sequence]
        return MSA(sequences, headers)

    @staticmethod
    def sum_of_pairs(msa, path_to_matrix) -> float:
        """
        Calculates alignment score using sum of pairs method

        :param msa:             MSA object of aligned sequences
        :param path_to_matrix:  path to the custom substitution matrix
        :return:                scoring of alignment
        """
        matrix = CustomMatrix(path_to_matrix)
        score = SumOfPairs(msa, matrix).compute()
        return score
