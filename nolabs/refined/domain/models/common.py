__all__ = ['Protein',
           'ProteinId',
           'ProteinContent',
           'ProteinName',
           'AminoAcidName',
           'AminoAcidId',
           'AminoAcid',
           'Job',
           'JobId',
           'JobName',
           'JobType']

import datetime
import os
from enum import Enum
from functools import cached_property
from typing import Union
from uuid import UUID

from mongoengine import EnumField, DateTimeField, EmbeddedDocument, FileField, Document, ReferenceField, CASCADE
from pydantic.dataclasses import dataclass

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.application.event_handlers.amino_acid_event_handlers import AminoAcidCreatedEvent
from nolabs.refined.application.event_handlers.protein_event_handlers import ProteinCreatedEvent
from nolabs.refined.domain.event_dispatcher import EventDispatcher
from nolabs.refined.domain.models.experiment import Experiment
from nolabs.refined.infrastructure.mongo_fields import ValueObjectUUIDField, ValueObjectStringField
from nolabs.seedwork.domain.entities import Entity
from nolabs.seedwork.domain.value_objects import ValueObject, ValueObjectBinary, ValueObjectString, ValueObjectUUID


@dataclass
class ExperimentId(ValueObjectUUID):
    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_experiment_id)

    def __eq__(self, other):
        if not isinstance(other, ProteinId):
            return False

        return self.value == other.value


@dataclass
class ProteinContent(ValueObjectBinary):
    def __post_init__(self):
        if not self.value:
            raise NoLabsException('Invalid protein content', ErrorCodes.invalid_protein_content)


@dataclass
class ProteinName(ValueObjectString):
    def __post_init__(self):
        if not self.value:
            raise NoLabsException('Protein name cannot be empty', ErrorCodes.invalid_protein_name)

        if len(self.value) < 1 or len(self.value) > 1000:
            raise NoLabsException('Length of protein name must be between 1 and 1000', ErrorCodes.invalid_protein_name)


@dataclass
class AminoAcidName(ValueObjectString):
    def __post_init__(self):
        if not self.value or not isinstance(self.value, str):
            raise NoLabsException('Amino acid name cannot be empty', ErrorCodes.invalid_aa_name)

        if '.fasta' in self.value:
            self.value = os.path.basename(self.value).replace('.fasta', '')

    @cached_property
    def fasta_name(self):
        return f'{self.value}.fasta'


@dataclass
class ExperimentName(ValueObjectString):
    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_experiment_name)

        if len(self.value) < 1 or len(self.value) > 100:
            NoLabsException.throw(ErrorCodes.invalid_experiment_name)


@dataclass
class ProteinId(ValueObjectUUID):
    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_protein_id)

    def __eq__(self, other):
        if not isinstance(other, ProteinId):
            return False

        return self.value == other.value


@dataclass
class AminoAcidId(ValueObjectUUID):
    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_aa_id)

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

    def __eq__(self, other):
        if not isinstance(other, ProteinId):
            return False

        return self.value == other.value


class JobName(ValueObjectString):
    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_job_name)

        if len(self.value) < 1 or len(self.value) > 100:
            NoLabsException.throw(ErrorCodes.invalid_job_name)


class JobType(str, Enum):
    LOCALISATION = 1
    FOLDING = 2


class AminoAcid(Document, Entity):
    id: AminoAcidId = ValueObjectUUIDField(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: AminoAcidName = ValueObjectStringField(required=True, unique=True)
    content = FileField(required=True)

    def __init__(self, id: AminoAcidId, experiment: Experiment, name: AminoAcidName, content: Union[bytes, str],
                 *args,
                 **values):
        if not id:
            raise NoLabsException('Amino acid id is invalid', ErrorCodes.invalid_aa_id)
        if not name:
            raise NoLabsException('Amino acid name is invalid', ErrorCodes.invalid_aa_name)
        if not content:
            raise NoLabsException('Amino acid content is invalid', ErrorCodes.invalid_aa_content)
        if not experiment:
            raise NoLabsException('Amino acid must be in experiment', ErrorCodes.invalid_experiment_id)

        super().__init__(id=id, experiment=experiment, name=name, content=content, *args, **values)

        EventDispatcher.raise_event(AminoAcidCreatedEvent(self))

    def set_name(self, value: AminoAcidName):
        if not value:
            NoLabsException.throw(ErrorCodes.invalid_aa_name)
        self.name = value


class Protein(Document, Entity):
    id: ProteinId = ValueObjectUUIDField(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: ProteinName = ValueObjectStringField(required=True)
    content = FileField(required=True)

    def __init__(
            self,
            id: ProteinId,
            experiment: Experiment,
            name: ProteinName,
            content: Union[bytes, str],
            *args,
            **kwargs
    ):
        if not id:
            raise NoLabsException('Protein id is invalid', ErrorCodes.invalid_protein_id)
        if not name:
            raise NoLabsException('Protein name is invalid', ErrorCodes.invalid_protein_name)
        if not experiment:
            raise NoLabsException('Protein must be in experiment', ErrorCodes.invalid_experiment_id)
        if not content:
            raise NoLabsException('Protein content is empty', ErrorCodes.invalid_protein_content)

        super().__init__(id=id, experiment=experiment, name=name, content=content, *args, **kwargs)

        EventDispatcher.raise_event(ProteinCreatedEvent(self))

    def set_name(self, value: ProteinName):
        if not value:
            NoLabsException.throw(ErrorCodes.invalid_protein_name)
        self.name = value


class Job(Document, Entity):
    id: JobId = ValueObjectUUIDField(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: JobName = ValueObjectStringField(required=True)
    job_type: JobType = EnumField(JobType, required=True)
    created_at: dataclass = DateTimeField(default=datetime.datetime.utcnow)
