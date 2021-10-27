"""
Unicode characters encoding/decoding
"""

import csv
import os
import re
import ntpath
from subprocess import run, PIPE, CompletedProcess
import sys


class AlignUtils:
    """
    Class of utility functions to handle replacement of Unicode characters to allow launching MAFFT
    in "text" mode to align Peptides with non-natural amino acids
    """

    def __init__(self):
        # class members
        self.__chars = []
        # dictionary for encoding
        self.__d_enc = {}
        # dictionary for decoding
        self.__d_dec = {}
        # list of natural amino acids
        self.__natAA = ["G", "A", "V", "L", "I", "M", "P", "F", "W", "S", "T", "N", "Q", "Y", "C", "K", "R", "H", "D", "E"]

    # Private methods:
    def __repl(self, val):
        """
        Replaces non-natural amino-acids to special symbols
        :param val: string, non-natural amino acid to be replaced by Hex code
        :return:
        """

        if len(self.__chars) == 0:
            exc_ascii = [0x3E, 0x3D, 0x3C, 0x2D, 0x20, 0x0d, 0x0a, 0x00]  # excluding ASCII chars
            all_chars = [chr(i) for i in range(256) if i not in exc_ascii]  # valid ASCII characters
            for x in all_chars:
                if x.upper() not in self.__natAA:
                    self.__chars.append(x)

        if val not in self.__d_enc:
            curr_pos = len(self.__d_enc)
            if curr_pos >= len(self.__chars):
                raise Exception("There are no more letters in mapping array")
            # assign ascii symbol one by one
            ch = self.__chars[curr_pos]
            self.__d_enc[val] = ch
            self.__d_dec[ch] = val
        result = self.__d_enc.get(val, '')  # if no match found, insert nothing
        return result

    # Instance methods:
    def encode_alignment_sequences(self, fasta_input_array, peptide_name, timestamp):
        """
        Replace non-natural AA
        NOTE: this method operates with ARRAY now, not with files
        :param fasta_input_array: array of FASTA input lines
        :param peptide_name: The Peptide name to be processed (PEPTIDE1, PEPTIDE2,...)
        :param timestamp: unique ID of this task
        :return:
        """

        output_file = os.path.join(peptide_name + "_" + str(timestamp) + "_mafft.txt")
        with open(output_file, "w", encoding='latin-1') as ofile:
            for name, seq in AlignUtils.read_fasta(fasta_input_array):  # calling generator f-ion
                # translated sequence according to sub matrix

                # Treat polymers of CHEM type as single monomer, encode the whole sequence
                if peptide_name.startswith("CHEM"):
                    monomer = re.match(r"\[.*\]", seq).group()
                    self.__repl(monomer)
                    seq = self.__d_enc[seq]
                else:
                    # findings = re.finditer(r'\[[^\]]*\]', seq)
                    findings = self.parse_naa_in_fasta(seq)
                    for s in findings:
                        # find all NAA replacements
                        self.__repl(s)

                    # replace NAA in original string
                    for key in self.__d_enc:
                        seq = seq.replace(key, self.__d_enc[key])
                ofile.write(name + "\n" + seq + "\n")

        return output_file

    def clear_symbols(self):
        """Clear private dictionaries of Unicode character codes"""
        self.__chars.clear()
        self.__d_enc.clear()
        self.__d_dec.clear()

    def decode_mafft(self, fasta_input_array):
        """
        Decodes MAFFT results
        :param fasta_input_array: FASTA input lines array
        :return:
        """
        fasta_output_array = []
        for name, seq in AlignUtils.read_fasta(fasta_input_array):
            translated = ""
            for b in seq:
                aa = self.__d_dec.get(b, b)
                translated += aa
            fasta_output_array.append(name + "\n" + translated + "\n")
        return fasta_output_array

    def create_substitution_matrix(self, matrix_file_path, monomers_map_file_path, peptide_name, timestamp):  # noqa: C901, pylint: disable=too-many-branches, too-many-statements
        """
        Create a custom substitution matrix
        NOTE: use this function during MAFFT calculation because it is required to provide a table of encoding for NAA
        Private variable __d_enc is set if nnAAs are present in the sequence, and not set if sequences contain only nAAs

        :param matrix_file_path: path to ROCS matrix file
        :param monomers_map_file_path: path to monomers map
        :param peptide_name: The Peptide name to be processed (PEPTIDE1, PEPTIDE2,...)
        :param timestamp: Unique marker for this API call
        :return:
        """

        # subst. matrix file
        subst_matrix_file = os.path.join(peptide_name + "_" + str(timestamp) + "_matrix.txt")

        # prepare list of used chars
        d_hex = {}
        for key in self.__d_enc:
            # remove brackets
            new_key = key.replace("[", "").replace("]", "")
            # hex(ord(self.__d_enc[key]))
            d_hex[new_key] = format(ord(self.__d_enc[key]), '#04x')
        for c in self.__natAA:
            d_hex[c] = hex(ord(c))
        # print(self.__d_hex)

        # read mapping from Monomers Map
        d_monomer_map = {}
        # list of non-natural AA which aren't found in Monomers Map
        not_found = []
        with open(monomers_map_file_path) as csv_file:
            reader = csv.reader(csv_file)
            header_skipped = False
            for row in reader:
                if not header_skipped:
                    header_skipped = True
                    continue  # skip header
                tokens = re.split(r'\t+', row[0])
                if len(tokens) == 0:
                    continue
                # Unicode field is key
                # array: Symbol field, position in mapping file
                symbol_field = tokens[0]
                unicode_field = AlignUtils.get_unicode_char(tokens[1])
                # aa = symbol_field.upper()
                if symbol_field in d_hex:
                    d_monomer_map[unicode_field] = [symbol_field, -1]
                    if len(d_hex) == len(d_monomer_map):
                        break  # all AA found and covered
        # print(d_monomer_map)
        if len(d_hex) != len(d_monomer_map):
            for key in d_hex:
                found = False
                for mm_key in d_monomer_map:
                    it = d_monomer_map[mm_key]
                    if it[0] == key:
                        found = True
                        break
                if not found:
                    not_found.append(key)
                    sys.stderr.write(key + " not found in Monomer Map\n")

        # reads ROCS file
        # indexing list of chars. The order of symbols are the same on vertical and horizontal
        with open(matrix_file_path, "rb") as mfp:
            # this row contains characters which we need for mapping with monomers_map
            # analyze the row to find coordinates
            first_row = mfp.readline()
            chars = first_row.split(b"\x20")
            # index = 0 has empty character. And this is good because rows also started with 1
            found_chars = 0
            # numbers of rows we need for quick access
            row_numbers = []
            for index in range(len(chars) - 1):
                current_char = chars[index].decode("utf-8")
                if current_char in d_monomer_map:
                    d_monomer_map[current_char][1] = index
                    row_numbers.append(index)
                    found_chars += 1
                    if found_chars == len(d_monomer_map):
                        break
                # print("index=" + str(index) + "; char=" + current_char)

        # check if we found all rows in ROCS
        if found_chars != len(d_monomer_map):
            for key in d_monomer_map:
                item = d_monomer_map[key]
                if item[1] == -1:
                    raise Exception("Can't find row with key=" + key + " in ROCS file")

        # lets find rows we need. We know their numbers. Other rows ignore.
        row_numbers.sort()
        last_row_number = row_numbers[len(row_numbers) - 1]
        rows_cache = {}
        # open file again for indexation
        with open(matrix_file_path, "rb") as mfp:
            for index, row in enumerate(mfp):
                if index in row_numbers:
                    arr = row.split(b"\x20")
                    ch = arr[0].decode("utf-8")
                    rows_cache[ch] = arr
                    if index == last_row_number:
                        break

        # prepare score matrix (char, char, score, comment)
        score_template = []
        for key_left in d_monomer_map:
            symbol_left = d_monomer_map[key_left][0]
            hex_left = d_hex[symbol_left]
            if key_left not in rows_cache:
                print("error")
            row_array = rows_cache[key_left]
            for key_right in d_monomer_map:
                # this index is 1-based but this is OK because the first position in a row is char
                dkr = d_monomer_map[key_right]
                symbol_right = dkr[0]
                hex_right = d_hex[symbol_right]
                column_index = dkr[1]
                score = row_array[column_index].decode("utf-8")
                score_template.append(str(hex_left) + " " + str(hex_right) + " " + str(
                    score) + "   # " + symbol_left + " x " + symbol_right)

        # add special rows for NAA which are not found in Monomer Map.
        if len(not_found) > 0:
            score_for__not_found = -10  # according to PEPSAR-59
            for naa in not_found:
                hex_left = d_hex[naa]
                for key_right in d_monomer_map:
                    # this index is 1-based but this is OK because the first position in a row is char
                    dkr = d_monomer_map[key_right]
                    symbol_right = dkr[0]
                    hex_right = d_hex[symbol_right]
                    score_template.append(str(hex_left) + " " + str(hex_right) + " " + str(
                        score_for__not_found) + "   # " + naa + " x " + symbol_right)

                # Add substitutions between not found monomers
                for naa2 in not_found:
                    score_for__not_found = -10
                    hex_right = d_hex[naa2]

                    # Set score to 10 in case of match
                    if hex_left == hex_right:
                        score_for__not_found = 10

                    score_template.append(str(hex_left) + " " + str(hex_right) + " " + str(
                        score_for__not_found) + "   # " + naa + " x " + naa2)

        with open(subst_matrix_file, "w", encoding='utf8') as ofl:
            for s in score_template:
                ofl.write(s + "\n")

        return subst_matrix_file

    # Static methods:
    @staticmethod
    def run_mafft_utility(*args, mafft_binary, matrix_file, gap_opening_penalty, gap_extension_penalty,
                          realign, realign_method="", mafft_options="") -> (str, str):
        """
        Run MAFFT utility
        :param args: list of input files for alignment with encoded characters
        :param mafft_binary: Absolute path to the MAFFT binary
        :param matrix_file: Matrix file. Use its name without path.
        :param gap_opening_penalty: a penalty for creating a gap on any length in an aligned sequence
        :param gap_extension_penalty: a penalty for extending a gap by one monomer. See https://en.wikipedia.org/wiki/Gap_penalty.
        :param realign: specifies, whether MAFFT should align sequences, or realign
        :param realign_method: MAFFT method used for realignment
        :param mafft_options: extra MAFFT options, copied unchanged to the MAFFT command string
        :return: MAFFT standard output and error as strings
        """

        head, tail = ntpath.split(matrix_file)
        matrix_file_fn_only = tail or ntpath.basename(head)
        matrix_part = " --textmatrix " + matrix_file_fn_only
        gaps = " --op " + str(gap_opening_penalty) + " --ep " + str(gap_extension_penalty)

        # If the shell locale is set to UTF-8, then MAFFT will not accept characters over 0x79,
        # and reports "tr: Illegal byte sequence"
        shell_locale_setting = "LC_CTYPE=C "
        mafft_cmd = shell_locale_setting + mafft_binary + matrix_part + gaps + " " + mafft_options + " --text "

        if realign:
            method = "--" + realign_method
            mafft_cmd += " " + method + " " + args[1] + " " + args[0]
        else:
            mafft_cmd += args[0]

        mafft_result: CompletedProcess = run(mafft_cmd, shell=True, check=True, stdout=PIPE, stderr=PIPE)

        return mafft_result.stdout.decode(encoding="latin-1"), mafft_result.stderr.decode(encoding="utf-8",
                                                                                          errors="ignore")

    @staticmethod
    def get_unicode_char(un):
        """
        Converts string to unicode character
        :param un: string representing the unicode character
        :return:
        """
        ln = len(un)
        if ln < 4:
            for _ in range(ln, 4):
                un = "0" + un
        return chr(int("0x" + un, 16))

    @staticmethod
    def read_fasta(fasta_input_array):
        """
        Read an input array of lines in FASTA format
        :param fasta_input_array: array of strings (FASTA lines)
        :return:
        """
        name, seq = None, []
        for item in fasta_input_array:
            tokens = item.split("\n")
            for line in tokens:
                line = line.rstrip()
                if len(line) == 0:
                    continue
                if line.startswith(">"):
                    if name:
                        yield name, "".join(seq)  # generator f-ion
                    name, seq = line, []
                else:
                    seq.append(line)
                    # seq.append(line.upper())  # converting seq to upper letters
        if name:
            yield name, "".join(seq)

    @staticmethod
    def parse_naa_in_fasta(sequence):
        """
        Parses FASTA sequence and extracts non-natural AAs, including monomers defined by SMILES

        :param sequence: FASTA sequence
        :return: list of non-natural AAs
        """
        naas = []

        # Track the number of opened and closed brackets, in case there are nested monomers, like SMILES
        opened = 0
        closed = 0

        naa = False
        monomer = ""

        for el in sequence:

            if el == "[":
                # If this is an opening bracket, then just increase opened variable and set naa to True
                if opened == 0:
                    naa = True
                    opened += 1
                else:
                    opened += 1
                    monomer += el

            elif el == "]":
                # If this is closing bracket, then append monomer to the list
                if opened == closed + 1:
                    naas.append("[" + monomer + "]")
                    opened = 0
                    closed = 0
                    monomer = ""
                    naa = False
                else:
                    monomer += el
                    closed += 1
            elif naa:
                monomer += el

        assert not naa, f"Non-valid notation of non-natural amino acids in the following sequence {sequence}"

        return naas
