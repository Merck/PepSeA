"""
API input/output models defined using Pydantic package
"""

from enum import Enum
from typing import List
from pydantic import BaseModel


# Models to be used in the API
# pylint: disable=missing-class-docstring
# Peptide model - '/align' API input
class Peptide(BaseModel):
    ID: str
    HELM: str

    class Config:
        schema_extra = {
            "example": {
                "ID": "L-000000001",
                "HELM": "PEPTIDE1{[ClAc].F.R.Y.L.Y.[Ahp].F.C.G.K.K.[NH2]}$PEPTIDE1,PEPTIDE1,1:R1-9:R3$$$V2.0"
            }
        }


# Peptide model - API output base
class AlignedPeptide(BaseModel):
    PolymerID: str
    AlignedSubpeptide: str
    HELM: str
    ID: str
    AlignedSeq: str


# API Status endpoint "/status" response model
class ApiStatus(BaseModel):
    Status: str

    class Config:
        schema_extra = {"example": {"Status": "OK"}}


# MAFFT version endpoint "/version" response model
class MafftVersion(BaseModel):
    MAFFT_version: str

    class Config:
        schema_extra = {"example": {"MAFFT_version": "v7.471 (2020/Jul/3)"}}


# MAFFT available alignment methods enumeration
class MafftMethods(str, Enum):
    auto = "mafft --auto"
    mafft = "mafft"
    linsi = "linsi"
    ginsi = "ginsi"
    einsi = "einsi"
    fftns = "fftns"
    fftnsi = "fftnsi"
    nwns = "nwns"
    nwnsi = "nwnsi"


# '/realign' endpoint input model
class RealignInput(BaseModel):
    aligned_sequences: List[AlignedPeptide]
    new_sequences: List[Peptide]


# MAFFT methods for realignment
class RealignMethods(str, Enum):
    add = "add"
    addfull = "addfull"
    addlong = "addlong"
    addfragments = "addfragments"
    addprofile = "addprofile"
