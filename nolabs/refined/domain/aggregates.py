from pydantic.dataclasses import dataclass

from nolabs.refined.domain.value_objects import AminoAcidName, ProteinName, AminoAcidContent, ProteinContent, ProteinId, \
    AminoAcidId
from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.seedwork.domain.aggregates import Aggregate


@dataclass
class AminoAcid(Aggregate):
    _id: AminoAcidId
    _name: AminoAcidName
    _content: AminoAcidContent

    def __init__(self, p_id: AminoAcidId, name: AminoAcidName, content: AminoAcidContent):
        if not p_id:
            NoLabsException.throw(ErrorCodes.invalid_aa_id)

        if not name:
            NoLabsException.throw(ErrorCodes.invalid_aa_name)

        if not content:
            NoLabsException.throw(ErrorCodes.invalid_aa_content)

        self._id = p_id
        self._name = name
        self._content = content

    @property
    def id(self) -> AminoAcidId:
        return self._id

    @property
    def name(self) -> AminoAcidName:
        return self._name

    def write_content(self, content: AminoAcidContent):
        if not content:
            NoLabsException.throw(ErrorCodes.invalid_aa_content)

        self._content = content

    def read_content(self) -> str:
        return self._content.read()


@dataclass
class Protein(Aggregate):
    _id: ProteinId
    _name: ProteinName
    _content: ProteinContent

    def __init__(self, p_id: ProteinId, name: ProteinName, content: ProteinContent):
        if not p_id:
            NoLabsException.throw(ErrorCodes.invalid_protein_id)

        if not name:
            NoLabsException.throw(ErrorCodes.invalid_protein_content)

        if not content:
            NoLabsException.throw(ErrorCodes.invalid_protein_content)

        self._id = p_id
        self._name = name
        self._content = content

    @property
    def id(self) -> ProteinId:
        return self._id

    @property
    def name(self) -> ProteinName:
        return self._name

    def write_content(self, content: ProteinContent):
        if not content:
            NoLabsException.throw(ErrorCodes.invalid_protein_content)

        self._content = content

    def read_content(self) -> str:
        return self._content.read()


@dataclass
class Job(Aggregate):

