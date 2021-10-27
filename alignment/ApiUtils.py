"""
API input/output manipulation utility functions
"""

import json

from alignment.PyHELM_simple import HelmObj


def extract_helm_from_json(input_json):
    """Extracts the HELM strings out of a JSON array.

    :param input_json: JSON array of Peptide objects
    :return: output_helm: string (extracted HELM strings separated by a newline)
    """
    output_helm = ""
    for peptide in input_json:
        if len(output_helm) > 0:
            output_helm += "\n"
        output_helm += peptide["HELM"]
    return output_helm


def json_output(helm_data, json_input):
    """
    Creates JSON output from the decoded and aligned HELM sequence and JSON input.
    Returns JSON array of aligned sequences.

    :param helm_data: Dictionary of the aligned subpeptides
    :param json_input: Original JSON input
    :return:
    """

    # If helm_data is empty, the function will return an empty json_list
    json_list = []
    if len(helm_data) > 0:  # pylint: disable=too-many-nested-blocks

        # The following cycle iterates over the helm_data dictionary, where keys are the names of subpeptides
        for key in helm_data:

            # peptides variable stores the sequences of subpeptides of the given key
            peptides = [x.strip("\n") for x in helm_data[key].split(f"> {key}\n")[1:]]
            index = 0

            """
            The following cycle iterates over the items of the original JSON input.

            At the beginning, it checks, whether the current entity contains the given subpeptide.
            If the subpeptide exists in current sequence, then "if" part is executed, which creates
            separate item for JSON output.


            """
            for entity in json_input:

                """
                Differentiate between aligned input (output of "/align" endpoint) and raw input.
                If we don't have already aligned input, we only need to extract polymer, which equals current key (subpeptide)
                If we have aligned input, we need to extract polymer, which equals key AND PolymerID
                """
                polymer_in_entity = entity.get('PolymerID', None)
                if not polymer_in_entity or polymer_in_entity == key:
                    sequence = extract_subpeptide(key, entity['HELM'])
                else:
                    sequence = None

                # If there is a subpeptide with a given key, we can add it to the output
                if sequence is not None:
                    output_json = "{"
                    output_json += '"PolymerID":"' + key + '", '
                    output_json += '"AlignedSubpeptide":"' + sequence + '", '
                    output_json += '"HELM":"' + entity['HELM'] + '", '
                    output_json += '"ID":"' + entity["ID"] + '", '
                    output_json += '"AlignedSeq":"' + peptides[index] + '", '

                    """
                    The following cycle iterates over the aligned subsequence, and the the help of the "extract_monomer"
                    function adds separate monomers to the JSON output in the following format:

                    "PEPTIDE1_1": "d1Nal",

                    """
                    if key.startswith("CHEM"):
                        output_json += f'"{key}_1":"' + sequence[1:-1] + '"' + ", "
                    else:
                        for ind, el in enumerate(extract_monomer(peptides[index])):
                            output_json += f'"{key}_{ind + 1}":"' + el + '"' + ", "
                    output_json = output_json[:-2]
                    output_json += "}"
                    index += 1

                    json_list.append(json.loads(output_json))

    return json_list


def extract_subpeptide(peptide_name, helm_string):
    """Extracts the sub-peptide sequence, given the HELM string and the name of the sub-peptide.

    :param peptide_name:
    :param helm_string:
    :return:
    """
    helm_obj = HelmObj()
    helm_obj.parse_helm(helm_string)

    # The following cycle goes through all the subpeptides of the HELM sequence.
    # If the name of the subpeptide equals the peptide_name, the function returns sequence of the subpeptide
    for poly in helm_obj.polymers:
        if poly.name == peptide_name:
            return poly.data
    return None


def extract_monomer(aligned_sequence):
    """
    Generator function. Yields the monomers of the aligned sequence. Every symbol is accounted for monomer,
    unless they are not in square brackets. Otherwise, symbols in square brackets are precessed as single monomer.

    Currently works with sequences, where monomers are not separated with dots, e.g.
        [ClAc]FRYLY[Ahp]FCGKK[NH2]

    If monomers will be separated with dots, it will be necessary to update this function.
    :param aligned_sequence:
    :return:
    """
    non_natural = False
    monomer = ''

    # These variables track the number of opened and and closed square brackets.
    # This is needed for correct parsing of sequences, which have SMILES as monomers
    opened = 0
    closed = 0

    for el in aligned_sequence:
        if monomer != "" and not non_natural:
            monomer = ""

        if el == "[" and opened == 0:
            non_natural = True
            opened += 1
        elif el == "[":
            opened += 1
            monomer += el
        elif not non_natural:
            yield el

        # Yield the monomer, only if the current symbol is closing bracket,
        # and number of opened brackets is exactly one more than the number of closed ones
        elif el == "]" and opened == closed + 1:
            non_natural = False
            opened = 0
            closed = 0
            yield monomer
        elif el == "]":
            monomer += el
            closed += 1
        else:
            monomer += el

    assert not non_natural, f"Non-valid notation of non-natural amino acids in the following sequence {aligned_sequence}"


def extract_aligned_sequences(aligned_array):
    """
    Transforms aligned sequences into FASTA format

    :param aligned_array: array of aligned sequences
    :return: string, aligned sequences in FASTA format
    """

    aligned_seqs = ''

    for element in aligned_array:
        if len(aligned_seqs) > 0:
            aligned_seqs += '\n'
        aligned_seqs += f"> {element['PolymerID']}\n"
        aligned_seqs += element['AlignedSeq']

    return aligned_seqs


def escape_double_quotes_in_input(input_data):
    """ Adds escape character in front of double quotes in HELM string

    :param input_data: JSON array of input data
    :return: JSON array
    """

    if "aligned_sequences" and "new_sequences" in input_data:

        for ind, element in enumerate(input_data["aligned_sequences"]):
            helm = element["HELM"].replace('"', '\\"')
            input_data["aligned_sequences"][ind]["HELM"] = helm

        for ind, element in enumerate(input_data["new_sequences"]):
            helm = element["HELM"].replace('"', '\\"')
            input_data["new_sequences"][ind]["HELM"] = helm

    else:

        for ind, element in enumerate(input_data):
            helm = element["HELM"].replace('"', '\\"')
            input_data[ind]["HELM"] = helm

    return input_data
