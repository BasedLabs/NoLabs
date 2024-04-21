__all__ = ['Protein',
           'ProteinId',
           'ProteinContent',
           'ProteinName',
           'AminoAcidName',
           'AminoAcidId',
           'AminoAcid',
           'AminoAcidContent',
           'Job',
           'JobId',
           'JobName',
           'JobState',
           'JobType']

import datetime
import os
from dataclasses import field
from enum import Enum
from functools import cached_property
from uuid import UUID

from pydantic.dataclasses import dataclass

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.seedwork.domain.aggregates import Aggregate
from nolabs.seedwork.domain.value_objects import ValueObject


@dataclass
class ProteinContent(ValueObject):
    value: bytes

    def __post_init__(self):
        if not self.value:
            raise NoLabsException('Invalid protein content', ErrorCodes.invalid_protein_content)


@dataclass
class ProteinName(ValueObject):
    _min_length = 1
    _max_length = 1000

    value: str

    def __post_init__(self):
        if not self.value:
            raise NoLabsException('Protein name cannot be empty', ErrorCodes.invalid_protein_name)

        if len(self.value) < self._min_length or len(self.value) > self._max_length:
            raise NoLabsException('Length of protein name must be between 1 and 1000', ErrorCodes.invalid_protein_name)


@dataclass
class AminoAcidContent(ValueObject):
    value: bytes

    def __post_init__(self):
        if not self.value or self.value[0] != b'>':
            raise NoLabsException('Invalid amino acid content', ErrorCodes.invalid_aa_content)

    def __str__(self):
        return self.value.decode('utf-8')


@dataclass
class AminoAcidName(ValueObject):
    value: str

    def __post_init__(self):
        if not self.value or not isinstance(self.value, str):
            raise NoLabsException('Amino acid name cannot be empty', ErrorCodes.invalid_aa_name)

        if '.fasta' in self.value:
            self.value = os.path.basename(self.value).replace('.fasta', '')

    @cached_property
    def fasta_name(self):
        return f'{self.value}.fasta'

@dataclass
class ExperimentName(ValueObject):
    value: str

    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_experiment_name)

        if len(self.value) < 1 or len(self.value) > 100:
            NoLabsException.throw(ErrorCodes.invalid_experiment_name)


@dataclass
class ProteinId(ValueObject):
    value: UUID

    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_protein_id)

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not isinstance(other, ProteinId):
            return False

        return self.value == other.value


@dataclass
class ExperimentId(ValueObject):
    value: UUID

    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_experiment_id)

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not isinstance(other, ProteinId):
            return False

        return self.value == other.value


@dataclass
class AminoAcidId(ValueObject):
    value: UUID

    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_aa_id)

    def __str__(self):
        return str(self.value)

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not isinstance(other, ProteinId):
            return False

        return self.value == other.value


@dataclass
class JobId(ValueObject):
    value: UUID

    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_job_id)

    def __str__(self):
        return str(self.value)

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not isinstance(other, ProteinId):
            return False

        return self.value == other.value


@dataclass
class JobName(ValueObject):
    value: str

    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_job_name)

        if len(self.value) < 1 or len(self.value) > 100:
            NoLabsException.throw(ErrorCodes.invalid_job_name)

    def __str__(self):
        return self.value


class JobState(str, Enum):
    CREATED = 'CREATED'
    STARTED = 'STARTED'
    STOPPED = 'STOPPED'
    ERROR = 'ERROR'
    FINISHED = 'FINISHED'


class JobType(str, Enum):
    LOCALISATION = 1
    FOLDING = 2



@dataclass
class AminoAcid(Aggregate):
    id: AminoAcidId
    name: AminoAcidName
    content: AminoAcidContent

    def __post_init__(self):
        if not self.id:
            raise NoLabsException('Amino acid id is invalid', ErrorCodes.invalid_aa_id)
        if self.name:
            raise NoLabsException('Amino acid name is invalid', ErrorCodes.invalid_aa_name)
        if not self.content:
            raise NoLabsException('Amino acid content is invalid', ErrorCodes.invalid_aa_content)

    def set_name(self, value: AminoAcidName):
        if not value:
            NoLabsException.throw(ErrorCodes.invalid_aa_name)
        self.name = value


@dataclass
class Protein(Aggregate):
    id: ProteinId
    name: ProteinName
    content: ProteinContent

    def __post_init__(self):
        if not self.id:
            raise NoLabsException('Protein id is invalid', ErrorCodes.invalid_protein_id)
        if not self.name:
            raise NoLabsException('Protein name is invalid', ErrorCodes.invalid_protein_name)
        if not self.content:
            raise NoLabsException

    def set_name(self, value: ProteinName):
        if not value:
            NoLabsException.throw(ErrorCodes.invalid_protein_name)
        self.name = value


@dataclass
class Job(Aggregate):
    id: JobId
    name: JobName
    job_type: JobType
    created_at: dataclass = field(default_factory=datetime.datetime.utcnow)
    state: JobState = field(default_factory=lambda: JobState.CREATED)

    def __post_init__(self):
        if not self.id:
            raise NoLabsException('Job id is invalid', ErrorCodes.invalid_job_id)
        if not self.name:
            raise NoLabsException('Job name is invalid', ErrorCodes.invalid_job_name)
        if not self.state:
            raise NoLabsException('Job state is invalid', ErrorCodes.invalid_job_state)
        if not self.job_type:
            raise NoLabsException('Job type is invalid', ErrorCodes.invalid_job_type)









