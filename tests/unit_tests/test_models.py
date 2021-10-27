import pytest

from alignment.models import Peptide, AlignedPeptide, ApiStatus, \
    MafftVersion, MafftMethods, RealignInput, RealignMethods


@pytest.fixture
def peptide_to_align():
    peptide = {"ID": "L-000000010",
               "HELM": "PEPTIDE1{H.[dS].Q.G.T.F.T.S.E.Y.S.K.Y.L.D.E.R.A.A.K.D.F.V.Q.W.L.L.N.K.[NH2]}|PEPTIDE2{[C18diacid].[gE].[AEEA].[AEEA]}$PEPTIDE1,PEPTIDE2,12:R3-4:R2|PEPTIDE1,PEPTIDE1,16:R3-20:R3$$$V2.0"
               }

    return peptide


@pytest.fixture
def aligned_peptide():
    peptide = {"PolymerID": "PEPTIDE1",
               "AlignedSubpeptide": "H.[dS].Q.G.T.F.T.S.E.Y.S.K.Y.L.D.E.R.A.A.K.D.F.V.Q.W.L.L.N.K.[NH2]",
               "HELM": "PEPTIDE1{H.[dS].Q.G.T.F.T.S.E.Y.S.K.Y.L.D.E.R.A.A.K.D.F.V.Q.W.L.L.N.K.[NH2]}|PEPTIDE2{[C18diacid].[gE].[AEEA].[AEEA]}$PEPTIDE1,PEPTIDE2,12:R3-4:R2|PEPTIDE1,PEPTIDE1,16:R3-20:R3$$$V2.0",
               "ID": "L-000000010",
               "AlignedSeq": "H--[dS]QGTFTSEYSKYLDERAAKDFVQWLLN----K[NH2]"
               }

    return peptide


"""
Below are test for Peptide() model
"""


def test_Peptide(peptide_to_align):
    """Test right initiation of Peptide() model"""

    peptide = Peptide(
        ID=peptide_to_align['ID'],
        HELM=peptide_to_align['HELM']
    )

    assert peptide.ID == peptide_to_align['ID']
    assert peptide.HELM == peptide_to_align['HELM']


def test_Peptide_validation_error():
    """Test raising of validation error on wrong initiation"""

    with pytest.raises(Exception) as exception:
        Peptide()

    assert exception.typename == "ValidationError"


"""
Below are test for AlignedPeptide() model
"""


def test_AlignedPeptide(aligned_peptide):
    """Test right initiation of AlignedPeptide() model"""

    output = AlignedPeptide(
        PolymerID=aligned_peptide["PolymerID"],
        AlignedSubpeptide=aligned_peptide['AlignedSubpeptide'],
        HELM=aligned_peptide['HELM'],
        ID=aligned_peptide['ID'],
        AlignedSeq=aligned_peptide['AlignedSubpeptide']
    )

    assert output.PolymerID == aligned_peptide["PolymerID"]
    assert output.AlignedSubpeptide == aligned_peptide['AlignedSubpeptide']
    assert output.HELM == aligned_peptide['HELM']
    assert output.ID == aligned_peptide['ID']
    assert output.AlignedSeq == aligned_peptide['AlignedSubpeptide']


def test_AlignedPeptide_validation_error():
    """Test raising of validation error on wrong initiation"""

    with pytest.raises(Exception) as exception:
        AlignedPeptide()

    assert exception.typename == "ValidationError"


"""
Below are test for ApiStatus() model
"""


def test_ApiStatus():

    value = "Ok"
    status = ApiStatus(Status=value)

    assert status.Status == value


def test_ApiStatus_validation_error():
    """Test raising of validation error on wrong initiation"""

    with pytest.raises(Exception) as exception:
        ApiStatus()

    assert exception.typename == "ValidationError"


"""
Below are test for MafftVersion() model
"""


def test_MafftVersion():

    value = "v7.471 (2020/Jul/3)"
    version = MafftVersion(MAFFT_version=value)

    assert version.MAFFT_version == value


def test_MafftVersion_validation_error():
    """Test raising of validation error on wrong initiation"""

    with pytest.raises(Exception) as exception:
        MafftVersion()

    assert exception.typename == "ValidationError"


"""
Below are test for MafftMethods() model
"""


def test_MafftMethods():
    methods = ["mafft --auto", "mafft", "linsi", "ginsi",
               "einsi", "fftns", "fftnsi", "nwns", "nwnsi"]

    for method in methods:
        mafft_method = MafftMethods(method)

        assert mafft_method == method


def test_MafftMethods_type_error():
    """Test raising of type error on wrong initiation"""

    with pytest.raises(Exception) as exception:
        MafftMethods()

    assert exception.typename == "TypeError"


"""
Below are test for RealignInput() model
"""


def test_RealignInput(peptide_to_align, aligned_peptide):

    realign_input = RealignInput(aligned_sequences=[aligned_peptide], new_sequences=[peptide_to_align])

    assert realign_input.aligned_sequences == [aligned_peptide]
    assert realign_input.new_sequences == [peptide_to_align]


def test_RealignInput_validation_error():
    """Test raising of validation error on wrong initiation"""

    with pytest.raises(Exception) as exception:
        RealignInput()

    assert exception.typename == "ValidationError"


"""
Below are test for RealignMethods() model
"""


def test_RealignMethods():

    methods = ["add", "addfull", "addlong", "addfragments", "addprofile"]

    for method in methods:
        realign_method = RealignMethods(method)
        assert realign_method == method


def test_RealignMethods_type_error():
    """Test raising of type error on wrong initiation"""

    with pytest.raises(Exception) as exception:
        RealignMethods()

    assert exception.typename == "TypeError"
