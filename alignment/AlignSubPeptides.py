"""
Main functions for Peptide MSA (Multiple Sequence Alignment)
"""

import os
import glob
from datetime import datetime
from typing import Dict

from alignment.AlignUtils import AlignUtils
from alignment.PyHELM_simple import HelmObj
from alignment.scoring.ScoringUtils import Scoring


def align_sub_peptides(*args, gap_opening_penalty, gap_extension_penalty, polymer_to_align,  # pylint: disable=too-many-locals
                       path_to_mafft, path_to_subst_matrix, path_to_monomer_table, mafft_options="",
                       realign_method="") -> (Dict[str, str], str, Dict[str, float]):
    """
    Align or realign (depends in input) HELM sequences with non-natural amino acids,
    separately for each sub-peptide (PEPTIDE1, PEPTIDE2, etc.).

    :param args: a list of HELM strings
    :param gap_opening_penalty: a penalty for creating a gap on any length in an aligned sequence
    :param gap_extension_penalty: a penalty for extending a gap by one monomer. See https://en.wikipedia.org/wiki/Gap_penalty.
    :param polymer_to_align: name of the subpeptide, used for alignment. If None, then all the subpeptides are aligned
    :param path_to_mafft: a path to the MAFFT binary
    :param path_to_subst_matrix: name of the file with the alignment substitution matrix
    :param path_to_monomer_table: name of the file wih a list of monomers, providing Unicode characters representing monomers in
                                  the substitution matrix
    :param mafft_options: extra MAFFT options, copied unchanged to the MAFFT call string
    :param realign_method: method for sequences realignment,
           options: ["add", "addfull", "addlong", "addfragments", "addprofile"]
    :return: 1) A dictionary of aligned sub-peptides: keys are sub-peptide names, and values are aligned sub-peptide sequences,
                one per line.
             2) The stderr output of MAFFT.
    """

    aligner = AlignUtils()

    # Unique marker of this task, used for creating temporary file names.
    timestamp = datetime.now().isoformat(timespec='milliseconds')

    # Checks the number of inputs. 1 input array means alignment, 2 means realignment
    realign = bool(len(args) == 2)

    if realign:
        # Create distinct file names for aligned sequences and new sequences
        timestamp_aligned = str(timestamp) + "_aligned"
        timestamp_new = str(timestamp) + "_new"

        # Transform aligned sequences into HELM strings
        aligned_seqs_helm = convert2helm(args[0])
        aligned_seqs_helm = "\n".join(aligned_seqs_helm)

        # Split HELM strings of aligned sequences into several sub-peptides
        aligned_seqs_helm_dict = split_sub_peptides(aligned_seqs_helm, polymer_to_align)

        # Split HELM strings of new sequences into several sub-peptides
        new_seqs_helm_dict = split_sub_peptides(args[1], polymer_to_align)

        # Get the list of subpeptides, which are present in both sets of sequences
        subpeptides = get_common_subpeptides(aligned_seqs_helm_dict.keys(), new_seqs_helm_dict.keys())
    else:
        # Split one HELM string to several sub-peptides
        helm_aligned_dict = split_sub_peptides(args[0], polymer_to_align)
        subpeptides = helm_aligned_dict.keys()

    # This will return error message, in case the specified PolymerID is not in the input
    if polymer_to_align is not None and polymer_to_align not in subpeptides:
        return None, None, None

    # Align all sub-peptides using MAFFT in a loop
    main_output = {}
    mafft_stderr = ""

    # Scores of alignments
    alignment_scores = {}

    for key in subpeptides:

        if realign:
            #  Convert aligned sequences into FASTA
            aligned_fasta_input_lines = []
            aligned_helm_input_lines = aligned_seqs_helm_dict[key].split("\n")
            for line in aligned_helm_input_lines:
                aligned_fasta_input_lines.append(helm2fasta(line))

            # Convert new sequences into FASTA
            new_fasta_input_lines = []
            new_helm_input_lines = new_seqs_helm_dict[key].split("\n")
            for line in new_helm_input_lines:
                new_fasta_input_lines.append(helm2fasta(line))

            # Create file with input data of aligned sequences for MAFFT utility
            # it contains FASTA rows with replaced nn-AA
            # here we store conversion info inside the aligner
            aligned_seqs_input_file = aligner.encode_alignment_sequences(aligned_fasta_input_lines,
                                                                         key,
                                                                         timestamp_aligned)

            """
            Create file with input data of new sequences for MAFFT utility
            It contains FASTA rows with replaced nn-AA
            Here we store conversion info inside the aligner
            The conversion info includes encoding of both aligned and new sequences
            This info is cleaned up at the end of the cycle
            """
            new_seqs_input_file = aligner.encode_alignment_sequences(new_fasta_input_lines,
                                                                     key,
                                                                     timestamp_new)
        else:
            #  Convert to FASTA
            fasta_input_lines = []
            helm_input_lines = helm_aligned_dict[key].split("\n")
            for line in helm_input_lines:
                fasta_input_lines.append(helm2fasta(line))

            # Skip, if there is only one sequence of the given sub-peptide
            if len(fasta_input_lines) == 1:
                main_output[key] = fasta_input_lines[0]
                alignment_scores[key] = None
                continue

            # Create file with input data for MAFFT utility
            # it contains FASTA rows with replaced nn-AA
            # here we store conversion info inside the aligner
            input_file_mafft = aligner.encode_alignment_sequences(fasta_input_lines, key, timestamp)

        # Calculate substitution matrix
        subst_matrix_file = aligner.create_substitution_matrix(path_to_subst_matrix,
                                                               path_to_monomer_table,
                                                               key,
                                                               timestamp)

        if realign:
            # Realign sequences in MAFFT
            mafft_stdout, mafft_stderr = aligner.run_mafft_utility(aligned_seqs_input_file, new_seqs_input_file,
                                                                   mafft_binary=path_to_mafft,
                                                                   matrix_file=subst_matrix_file, gap_opening_penalty=gap_opening_penalty,
                                                                   gap_extension_penalty=gap_extension_penalty, realign=realign,
                                                                   realign_method=realign_method, mafft_options=mafft_options)
        else:
            # Run MAFFT utility
            mafft_stdout, mafft_stderr = aligner.run_mafft_utility(input_file_mafft, mafft_binary=path_to_mafft,
                                                                   matrix_file=subst_matrix_file,
                                                                   gap_opening_penalty=gap_opening_penalty,
                                                                   gap_extension_penalty=gap_extension_penalty,
                                                                   realign=realign, mafft_options=mafft_options)

        # Read aligned sequences in FASTA from the full MAFFT output
        fasta_mafft_output_array = get_aligned_sequences(mafft_stdout)

        # Calculate alignment score
        alignment_scores[key] = get_alignment_score(mafft_stdout, subst_matrix_file)

        # Decode MAFFT output and return the real names for replaced NAA
        fasta_mafft_decoded_array = aligner.decode_mafft(fasta_mafft_output_array)

        # Convert array to string
        main_output[key] = "".join(fasta_mafft_decoded_array)

        # Clean data about encoded sequences from the aligner
        aligner.clear_symbols()

    # Clean-up temporary files for this task
    cleanup(timestamp)

    return main_output, mafft_stderr, alignment_scores


def split_sub_peptides(helm_input, polymer_to_align):
    """
    Split HELM string into individual sub-peptides (PEPTIDE1, PEPTIDE2...)

    :param helm_input: HELM string
    :param polymer_to_align: name of the subpeptide, used for alignment. If None, then all the subpeptides are aligned
    :return: helm_lines: dictionary [PEPTIDE_NAME, HELM_STR] (instead peptide may be any other key)
    """
    input_lines = helm_input.split("\n")
    helm_obj = HelmObj()
    helm_lines = {}
    for line in input_lines:
        helm_obj.parse_helm(line)
        for polymer in helm_obj.polymers:
            if not polymer_to_align or polymer.name == polymer_to_align:
                helm_lines = extract_sub_peptide(polymer.name, helm_lines, polymer.data)

    return helm_lines


def extract_sub_peptide(polymer_name, helm_dict, polymer_data):
    """
    Extracts the specified subpeptide and adds it to the dictionary of HELM strings

    :param polymer_name: name of the polymer to extract
    :param helm_dict: dictionary of the HELM strings, where keys are the names of the polymers
    :param polymer_data: HELM sequence of the given subpeptide
    :return: updated dictionary of the HELM strings
    """
    if polymer_name not in helm_dict:
        helm_dict[polymer_name] = ""
    v = helm_dict[polymer_name]
    if len(v) > 0:
        v += "\n"
    v += polymer_name + "{" + polymer_data + "}$$$$"
    helm_dict[polymer_name] = v

    return helm_dict


def helm2fasta(helm):
    """
    Converts HELM string into an extended FASTA format.
    Only peptide polymers are converted, other polymers are ignored.
    Polymer names become headers in the FASTA output.
    Any residue names longer than one character are converted into "(name)" in the FASTA output.

    :param helm: HELM string
    :return: fasta: FASTA string
    """
    helm_obj = HelmObj()
    helm_obj.parse_helm(helm)

    fasta = ''
    for polymer in helm_obj.polymers:
        fasta += "> %s\n" % polymer.name

        # Read polymers of CHEM type as single monomer
        if polymer.type == "CHEM":
            fasta += polymer.data
        else:
            for monomer in helm_obj.get_monomers_from_polymer(helm, polymer.name):
                if len(monomer) == 1:
                    fasta += monomer
                else:
                    if "[" in monomer:
                        fasta += monomer
                    else:
                        fasta += "[%s]" % monomer
        fasta += '\n'

    return fasta


def convert2helm(fasta):
    """ Converts data from FASTA format to HELM string.
        FASTA headers become polymer type.
        Input is content of FASTA file in form of list
        Outputs a list of HELM strings
    """
    helm_list = []
    helm = ''

    fasta_list = fasta.split('\n')

    assert len(fasta_list) > 1 and fasta.startswith(">"), f"Non-valid FASTA notation:  {fasta}"

    for line in fasta_list:
        if line.startswith('>'):
            helm += line.strip('> \n')
        else:
            line = line.strip('\n')
            line = fasta2helm(line)
            helm += "{" + f"{line}" + "}$$$$"
            helm_list.append(helm)
            helm = ''
    return helm_list


def fasta2helm(fasta_seq):
    """Converts multi-line FASTA format into HELM strings"""
    helm = ''
    nnAA = False
    amino_symbol = ''
    for el in fasta_seq:
        if len(helm) > 0 and not nnAA:
            helm += "."

        if el == "[":
            amino_symbol += el
            nnAA = True
        elif el == "]":
            amino_symbol += el
            helm += amino_symbol
            nnAA = False
            amino_symbol = ''
        elif nnAA:
            amino_symbol += el
        else:
            helm += el

    # Check that all the nnAAs were processed correctly - all the squared brackets in the input sequence are paired
    assert not nnAA, f"Non-valid notation of non-natural amino acids in the following sequence {fasta_seq}"

    return helm


def get_common_subpeptides(list1, list2):
    """Get intersection of two sub-peptide chain lists (chain names)"""
    common_list = [sub_peptide for sub_peptide in list1 if sub_peptide in list2]
    return common_list


def get_aligned_sequences(mafft_output):
    """
    Parse aligned FASTA sequences from MAFFT output.
    :param mafft_output: MAFFT program output in FASTA format
    :return: Array of the aligned sequences in FASTA format
    """
    # mafft_lines = Array of FASTA lines

    mafft_lines = mafft_output.splitlines()
    headers = [line for line in mafft_lines if line.startswith(">")]
    sequences = ["".join(sequence.split("\n")[1:]) for sequence in mafft_output.split(">") if sequence]

    sequence_part = [i + "\n" + j + "\n" for i, j in zip(headers, sequences)]

    return sequence_part


def cleanup(timestamp):
    """
    Delete all temporary files
    :param timestamp:
    :return:
    """
    files = glob.glob(f'*{timestamp}*')
    for file in files:
        os.remove(file)


def get_alignment_score(mafft_out, path_to_subst_matrix):
    """
    Calculate the score of MSA

    :param mafft_out:               output from MAFFT program
    :param path_to_subst_matrix:    path to the substitution matrix
    :return:                        scoring of alignment
    """

    # Transform aligned sequences to MSA object
    # It is necessary as methods of pyMSA package work with these objects as input, rather than raw sequences
    msa = Scoring.mafft_output_to_msa(mafft_out)

    """
    In some cases there is only one polymer of a given type in the input data.
    Then, there is no possibility to evaluate the alignment.

    So algorithm will only calculate score for 2 or more aligned sequences, otherwise return None

    """

    if msa.number_of_sequences > 1:
        score = Scoring.sum_of_pairs(msa, path_to_subst_matrix)

        # Normalize the score by the number of sequences and the length of the alignment
        normalized_score = score / msa.number_of_sequences / msa.__len__()

        return normalized_score
    return None
