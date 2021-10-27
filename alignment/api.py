"""
Main API module - FastAPI endpoints
"""
# pylint: disable=line-too-long

import os
import subprocess
from typing import List, Optional

from fastapi import FastAPI, Query
from fastapi.encoders import jsonable_encoder
from fastapi.responses import RedirectResponse, PlainTextResponse, JSONResponse
from fastapi.middleware.httpsredirect import HTTPSRedirectMiddleware
from fastapi.middleware.cors import CORSMiddleware

from alignment.AlignSubPeptides import align_sub_peptides
from alignment.ApiUtils import extract_helm_from_json, json_output, extract_aligned_sequences, \
    escape_double_quotes_in_input
from alignment.models import Peptide, MafftMethods, ApiStatus, MafftVersion, RealignInput, \
    RealignMethods


# Environment variables
rocs_path = os.getenv("ROCS_FILE", "alignment/ROCS")
monomers_map_file = os.getenv("MONOMERS_MAP_FILE", "alignment/monomers_map.txt")
mafft_dir = os.getenv("MAFFT_DIR", "alignment/mafft/bin/")
mafft_man_file = os.getenv("MAFFT_MAN_FILE", "alignment/mafft_manpage.txt")
add_https_middleware = os.getenv("ADD_HTTPS_MIDDLEWARE", True)

# OpenAPI tags
tags_metadata = [{"name": "alignment", "description": "Multiple sequence alignment."}]

# List of CORS origins
# For now we use "*" wildcard - CORS enabled for all origins, but this disables use of credentials
# https://fastapi.tiangolo.com/tutorial/cors/
origins = [
    "*"
]


# Initialize FastAPI object
api = FastAPI(
    title="Peptide SAR",
    description='''This is an API for multiple sequence alignment based on the MAFFT program:
                <a href="https://mafft.cbrc.jp/alignment/software">https://mafft.cbrc.jp/alignment/software</a>.<br><br>
                It can also align peptide sequences containing non-natural amino acids.''',
    version="1.0.0",
    openapi_tags=tags_metadata
)


# Add middleware for HTTP to HTTPS redirection if deploying to Cloud Foundry
if add_https_middleware:
    api.add_middleware(HTTPSRedirectMiddleware)

# Add CORS middleware to enable requests coming from Spotfire
api.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)


# API Base URL: Redirect to Swagger documentation page
@api.get('/', include_in_schema=False)
def index():
    """API Base URL - redirect to OpenAPI documentation"""
    response = RedirectResponse(url='/docs')
    return response


@api.post('/align', tags=["alignment"], summary="Multiple sequence alignment.")
def align_sequences(peptides: List[Peptide],
                    method: MafftMethods = Query(default=MafftMethods.ginsi,
                                                 title="MAFFT alignment method",
                                                 description='''MAFFT alignment method used.<br>
                                                             See "/man" endpoint or MAFFT manual for more information.'''),
                    gap_open: float = Query(..., title="Gap opening penalty.",
                                            description="Penalty to open a gap of any length used in the global alignment.<br>"
                                                        "MAFFT default value: 1.53"),
                    gap_extend: float = Query(..., title="Gap extension penalty.",
                                              description="Penalty to extend an existing gap in the global alignment.<br>"
                                                          "MAFFT default value: 0.0"),
                    polymer_id: Optional[str] = Query(None, title="Polymer selected for alignment.",  # pylint: disable=unsubscriptable-object
                                                      description='Specify which chain will be used for alignment, example: "**PEPTIDE1**".'
                                                                  '<br>Default is **None** (All peptide chains found will be aligned).'),
                    mafft_options: str = Query(default="", title="Additional MAFFT arguments.",
                                               description="Extra arguments that can be passed to the MAFFT command line")):
    """Align HELM sequences of peptides with non-natural amino acids.

    \f
    :param peptides: JSON array of Peptide objects.
    :param method: MAFFT alignment method,
           enumeration: ["mafft --auto", "mafft", "linsi", "ginsi", "einsi", "fftnsi", "fftns", "nwns", "nwnsi"]
    For more information see "man mafft"
    :param gap_open: Gap opening penalty.
    :param gap_extend: Gap extension penalty.
    :param polymer_id: What polymer will be used for alignment, Example: PEPTIDE1.
    :param mafft_options: Extra MAFFT arguments.
    :return:

    Input example for body:

        [ ...,
        { 'ID': 'L-006424318-000P',
          'HELM': 'PEPTIDE1{[ClAc].F.F.Y.G.P.W.S.E.W.Y.F.Q.A.C.G.[NH2]}|PEPTIDE2{A.C.A.K.C.A}$PEPTIDE1,PEPTIDE2,2:R3-2:R2$$$V2.0'
          },
          ..., ]

    Input specifications:

    Key:                Value:
    ID                  ID of compound, any string, not null
    HELM                HELM string, not null

    Output example:

      { 'Alignment' :
                  [ ...,
                  { 'PolymerID': 'PEPTIDE1'
                    'AlignedSubpeptide': '[ClAc].F.F.Y.G.P.W.S.E.W.Y.F.Q.A.C.G.[NH2]',
                    'HELM': 'PEPTIDE1{[ClAc].F.F.Y.G.P.W.S.E.W.Y.F.Q.A.C.G.[NH2]}|PEPTIDE2{A.C.A.K.C.A}$PEPTIDE1,PEPTIDE2,2:R3-2:R2$$$V2.0',
                    'ID': 'L-006424318-000P',
                    'AlignedSeq: '[ClAc]F-FYGPWSEWYFQAC-G--[NH2]',
                    'PEPTIDE1_1': 'ClAc',
                    'PEPTIDE1_2': 'F',
                    'PEPTIDE1_3': '-',
                    'PEPTIDE1_4': 'F',
                    'PEPTIDE1_5': 'Y',
                    'PEPTIDE1_6': 'G',
                    'PEPTIDE1_7': 'P',
                    'PEPTIDE1_8': 'W',
                    'PEPTIDE1_9': 'S',
                    'PEPTIDE1_10': 'E',
                    'PEPTIDE1_11': 'W',
                    'PEPTIDE1_12': 'Y',
                    'PEPTIDE1_13': 'F',
                    'PEPTIDE1_14': 'Q',
                    'PEPTIDE1_15': 'A',
                    'PEPTIDE1_16': 'C',
                    'PEPTIDE1_17': '-',
                    'PEPTIDE1_18': 'G',
                    'PEPTIDE1_19': '-',
                    'PEPTIDE1_20': '-',
                    'PEPTIDE1_21': 'NH2'}
                      ..., ],

        'AlignmentScore' : 123987 }

    Output specifications:

    Alignment:
      Parameter           Description
      PolymerID           ID of the aligned polymer
      AlignedSubpeptide   Original sequence of the aligned subpeptide
      ID                  Input ID string
      HELM                Input HELM string
      AlignedSeq          Aligned sequence
      POLYMERID_1         Aligned monomer at position 1
      POLYMERID_i         Aligned monomer at position i

    AlignmentScore:
      Score of alignment normalized by the number of sequences and alignment's length
    """

    # Select which MAFFT alignment method to use (default is: "ginsi")
    mafft_binary = mafft_dir + method

    # This is necessary because input is of type Pydantic BaseModel object, which can be converted to dict
    # https://fastapi.tiangolo.com/tutorial/encoder/
    peptides_encoded = escape_double_quotes_in_input(jsonable_encoder(peptides))
    helm_input = extract_helm_from_json(peptides_encoded)

    if polymer_id:
        polymer_id = polymer_id.upper()

    align_output, _, alignment_score = align_sub_peptides(helm_input,
                                                          gap_opening_penalty=gap_open,
                                                          gap_extension_penalty=gap_extend,
                                                          polymer_to_align=polymer_id,
                                                          mafft_options=mafft_options,
                                                          path_to_mafft=mafft_binary,
                                                          path_to_subst_matrix=rocs_path,
                                                          path_to_monomer_table=monomers_map_file)
    if align_output is None:
        message = {"Message": f"Alignment failure. "
                              f"There are no sub-peptides with ID: {polymer_id} in the given input data."}
        return JSONResponse(status_code=460, content=message)

    output = json_output(align_output, peptides_encoded)

    return {"Alignment": output, "AlignmentScore": alignment_score}


@api.post('/realign', tags=["alignment"], summary="Multiple sequence realignment.")
def realign_sequences(input_json: RealignInput,
                      method: MafftMethods = Query(default=MafftMethods.ginsi,
                                                   title="MAFFT alignment method",
                                                   description='''MAFFT alignment method used.<br>
                                                               See "/man" endpoint or MAFFT manual for more information.'''),
                      realign_method: RealignMethods = Query(default=RealignMethods.add,
                                                             title="MAFFT realignment method",
                                                             description="Method for adding new sequences into existing alignment"),
                      gap_open: float = Query(..., title="Gap opening penalty.",
                                              description="Penalty to open a gap of any length used in the global alignment.<br>"
                                                          "MAFFT default value: 1.53"),
                      gap_extend: float = Query(..., title="Gap extension penalty.",
                                                description="Penalty to extend an existing gap in the global alignment.<br>"
                                                            "MAFFT default value: 0.0"),
                      polymer_id: Optional[str] = Query(None, title="Polymer selected for alignment.",  # pylint: disable=unsubscriptable-object
                                                        description='Specify which chain will be used for alignment,'
                                                                    'example: "**PEPTIDE1**".'
                                                                    '<br>Default is **None** (All peptide chains found will be aligned).'),
                      mafft_options: str = Query(default="", title="Additional MAFFT arguments.",
                                                 description="Extra arguments that can be passed to the MAFFT command line")
                      ):
    """
    Realign HELM sequences of peptides with non-natural amino acids.

    \f
    :param input_json: JSON,which includes "aligned_sequences" and "new_sequences" entries. Both are JSON arrays.
    :param method: MAFFT alignment method,
           enumeration: ["mafft --auto", "mafft", "linsi", "ginsi", "einsi", "fftnsi", "fftns", "nwns", "nwnsi"]
    For more information see "man mafft"
    :param realign_method: method for sequences realignment, options: ["add", "addfull", "addlong", "addfragments", "addprofile"]
    :param gap_open: Gap opening penalty.
    :param gap_extend: Gap extension penalty.
    :param polymer_id: What polymer will be used for alignment, Example: PEPTIDE1.
    :param mafft_options: Extra MAFFT arguments.
    :return:

    Input example for body:

        { 'aligned_sequences':
              [ ...,
              { 'PolymerID': 'PEPTIDE1'
                'AlignedSubpeptide': '[ClAc].F.F.Y.G.P.W.S.E.W.Y.F.Q.A.C.G.[NH2]',
                'HELM': 'PEPTIDE1{[ClAc].F.F.Y.G.P.W.S.E.W.Y.F.Q.A.C.G.[NH2]}|PEPTIDE2{A.C.A.K.C.A}$PEPTIDE1,PEPTIDE2,2:R3-2:R2$$$V2.0',
                'ID': 'L-006424318-000P',
                'AlignedSeq: '[ClAc]F-FYGPWSEWYFQAC-G--[NH2]',
                'PEPTIDE1_1': 'ClAc',
                'PEPTIDE1_2': 'F',
                'PEPTIDE1_3': '-',
                'PEPTIDE1_4': 'F'
                },
                ... ],

        'new_sequences' :
             [  ...,
                { 'ID': 'L-006424318-000P',
                'HELM': 'PEPTIDE1{[ClAc].F.F.Y.G.P.W.S.E.W.Y.F.Q.A.C.G.[NH2]}|PEPTIDE2{A.C.A.K.C.A}$PEPTIDE1,PEPTIDE2,2:R3-2:R2$$$V2.0'
                },
                ... ]
        }
    Input specifications:

    {aligned_sequences:
      Key:                Value:
      PolymerID           ID of the aligned polymer
      AlignedSubpeptide   Original sequence of the aligned subpeptide
      ID                  Input ID string
      HELM                Input HELM string
      AlignedSeq          Aligned sequence
      POLYMERID_1         Aligned monomer at position 1
      POLYMERID_i         Aligned monomer at position i

    new_sequences:
      Key:      Value:
        ID      ID of compound, any string, not null
        HELM    HELM string, not null
    }

    Output is in the same format, as in /align endpoint

    """

    # Select which MAFFT alignment method to use (default is: "ginsi")
    mafft_binary = mafft_dir + method

    if polymer_id:
        polymer_id = polymer_id.upper()

    decoded_input = escape_double_quotes_in_input(jsonable_encoder(input_json))

    # Extract aligned sequences from JSON
    aligned_sequences = extract_aligned_sequences(decoded_input['aligned_sequences'])

    # Extract new sequences
    new_helm_sequences = extract_helm_from_json(decoded_input['new_sequences'])

    realign_output, _, realign_score = align_sub_peptides(aligned_sequences,
                                                          new_helm_sequences,
                                                          gap_opening_penalty=gap_open,
                                                          gap_extension_penalty=gap_extend,
                                                          polymer_to_align=polymer_id,
                                                          path_to_mafft=mafft_binary,
                                                          path_to_subst_matrix=rocs_path,
                                                          path_to_monomer_table=monomers_map_file,
                                                          mafft_options=mafft_options,
                                                          realign_method=realign_method)

    if realign_output is None:
        message = {"Message": f"Alignment failure. "
                              f"There are no sub-peptides with ID: {polymer_id} in the given input data."}
        return JSONResponse(status_code=460, content=message)

    full_input_list = [*decoded_input['aligned_sequences'], *decoded_input['new_sequences']]

    output = json_output(realign_output, full_input_list)

    return {"Alignment": output, "AlignmentScore": realign_score}


@api.get('/version', response_model=MafftVersion, tags=["alignment"])
def mafft_version():
    """Print MAFFT version and exit."""
    result = subprocess.run([f'{mafft_dir + "mafft"}', '--version'], check=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    return {'MAFFT_version': result.stdout.decode().rstrip()}


@api.get('/status', response_model=ApiStatus, tags=["alignment"])
def api_status():
    """API health check endpoint for monitoring and alerting."""
    return {'Status': 'OK'}


@api.get('/man', response_class=PlainTextResponse, tags=["alignment"])
def mafft_manual():
    """Returns the whole MAFFT manual page as a plain text."""

    mafft_manual_output = ''

    with open(mafft_man_file, 'r', encoding="utf-8") as f:
        for line in f.readlines():
            mafft_manual_output += line

    return mafft_manual_output


@api.get('/help', response_class=PlainTextResponse, tags=["alignment"])
def mafft_help():
    """Returns the output of "mafft --help" command as a plain text."""

    mafft_help_output = subprocess.run(f'! {mafft_dir + "mafft"} --help', shell=True, check=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return mafft_help_output.stdout.decode()
