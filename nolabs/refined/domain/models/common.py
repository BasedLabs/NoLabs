__all__ = ['Protein',
           'ProteinId',
           'ProteinName',
           'AminoAcidName',
           'AminoAcidId',
           'AminoAcid',
           'Experiment',
           'Job',
           'JobId',
           'JobName',
           'AminoAcidCreatedEvent',
           'LocalisationProbability',
           'ProteinCreatedEvent']

import datetime
import os
import uuid
from functools import cached_property
from typing import Union
from uuid import UUID

from mongoengine import DateTimeField, Document, ReferenceField, CASCADE, EmbeddedDocument, \
    FloatField, EmbeddedDocumentField, BinaryField, UUIDField
from pydantic import model_validator
from pydantic.dataclasses import dataclass
from typing_extensions import Self

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.event_dispatcher import EventDispatcher
from nolabs.refined.infrastructure.mongo_fields import ValueObjectStringField
from nolabs.seedwork.domain.entities import Entity
from nolabs.seedwork.domain.events import DomainEvent
from nolabs.seedwork.domain.value_objects import ValueObject, ValueObjectString, ValueObjectUUID


@dataclass
class ExperimentId(ValueObjectUUID):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)
        return self

    def __eq__(self, other):
        if not isinstance(other, ExperimentId):
            return False

        return self.value == other.value


@dataclass
class ExperimentName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_experiment_name)

        if len(self.value) < 1 or len(self.value) > 1000:
            raise NoLabsException(ErrorCodes.invalid_experiment_name)
        return self


class Experiment(Document, Entity):
    id: UUID = UUIDField(primary_key=True, required=True)
    name: ExperimentName = ValueObjectStringField(required=True, factory=ExperimentName)
    created_at: datetime.datetime = DateTimeField(default=datetime.datetime.utcnow)

    def __init__(self, id: ExperimentId, name: ExperimentName, created_at: datetime.datetime | None = None, *args,
                 **kwargs):
        if not id:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)

        if not name:
            raise NoLabsException(ErrorCodes.invalid_experiment_name)

        created_at = created_at if created_at else datetime.datetime.now(tz=datetime.timezone.utc)

        super().__init__(id=id.value if isinstance(id, ExperimentId) else id, name=name, created_at=created_at, *args, **kwargs)

    @property
    def iid(self) -> ExperimentId:
        return ExperimentId(self.id)

    def set_name(self, name: ExperimentName):
        if not name:
            raise NoLabsException(ErrorCodes.invalid_experiment_name)

        self.name = name


@dataclass
class ProteinName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_protein_name)

        if len(self.value) < 1 or len(self.value) > 1000:
            raise NoLabsException(ErrorCodes.invalid_protein_name)
        return self


@dataclass
class AminoAcidName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value or not isinstance(self.value, str):
            raise NoLabsException(ErrorCodes.invalid_aa_name)

        if '.fasta' in self.value:
            self.value = os.path.basename(self.value).replace('.fasta', '')
        return self

    @cached_property
    def fasta_name(self):
        return f'{self.value}.fasta'


@dataclass
class ProteinId(ValueObjectUUID):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_protein_id)
        return self


@dataclass
class AminoAcidId(ValueObjectUUID):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_aa_id)
        return self

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not isinstance(other, AminoAcidId):
            return False

        return self.value == other.value


@dataclass
class JobId(ValueObjectUUID):
    value: UUID

    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_job_id)
        return self


@dataclass
class JobName(ValueObjectString):
    @model_validator(mode='after')
    def post_root(self) -> Self:
        if not self.value:
            raise NoLabsException(ErrorCodes.invalid_job_name)

        if len(self.value) < 1 or len(self.value) > 100:
            raise NoLabsException(ErrorCodes.invalid_job_name)
        return self


class LocalisationProbability(EmbeddedDocument, ValueObject):
    cytosolic: float = FloatField(required=True)
    mitochondrial: float = FloatField(required=True)
    nuclear: float = FloatField(required=True)
    other: float = FloatField(required=True)
    extracellular: float = FloatField(required=True)

    def __init__(self, cytosolic: float, mitochondrial: float,
                 nuclear: float, other: float, extracellular: float, *args, **kwargs):
        values = [
            cytosolic,
            mitochondrial,
            nuclear,
            other,
            extracellular
        ]

        for value in values:
            if not value:
                raise NoLabsException(ErrorCodes.invalid_localisation_probability)
            if value < 0 or value > 1.0:
                raise NoLabsException(ErrorCodes.invalid_localisation_probability)

        super().__init__(
            cytosolic=cytosolic,
            mitochondrial=mitochondrial,
            nuclear=nuclear,
            other=other,
            extracellular=extracellular,
            *args,
            **kwargs)


class AminoAcid(Document, Entity):
    id: UUID = UUIDField(primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: AminoAcidName = ValueObjectStringField(required=True, factory=AminoAcidName)
    content = BinaryField(required=True)

    localisation: LocalisationProbability = EmbeddedDocumentField(LocalisationProbability, required=False)

    def __init__(self, id: AminoAcidId, experiment: Experiment, name: AminoAcidName, content: Union[bytes, str], *args,
                 **kwargs):
        if not id:
            raise NoLabsException(ErrorCodes.invalid_aa_id)
        if not name:
            raise NoLabsException(ErrorCodes.invalid_aa_name)
        if not content:
            raise NoLabsException(ErrorCodes.invalid_aa_content)
        if not experiment:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)

        if isinstance(content, str):
            content = content.encode('utf-8')

        super().__init__(id=id.value if isinstance(id, AminoAcidId) else id, experiment=experiment, name=name, content=content, *args, **kwargs)

        EventDispatcher.raise_event(AminoAcidCreatedEvent(self))

    @property
    def iid(self) -> AminoAcidId:
        return AminoAcidId(self.id)

    def set_name(self, name: AminoAcidName):
        if not name:
            raise NoLabsException(ErrorCodes.invalid_aa_name)
        self.name = name

    def set_content(self, content: Union[bytes, str]):
        if not content:
            raise NoLabsException(ErrorCodes.invalid_aa_content)

        if isinstance(content, str):
            content = content.encode('utf-8')

        self.content = content

    def get_content(self) -> str | None:
        if self.content:
            return self.content.decode('utf-8')

        return None

    def set_localisation_probability(self, localisation: LocalisationProbability):
        if not localisation:
            raise NoLabsException(ErrorCodes.invalid_localisation_probability)

        self.localisation = localisation

    @classmethod
    def create(cls, experiment: Experiment, name: AminoAcidName, content: Union[bytes, str]):
        amino_acids = cls.objects.filter(name=name.value, experiment=experiment)
        if amino_acids:
            amino_acid: AminoAcid = amino_acids[0]
            amino_acid.set_content(content)
            amino_acid.set_name(name)
            return amino_acid

        return AminoAcid(
            AminoAcidId(uuid.uuid4()),
            experiment=experiment,
            name=name,
            content=content
        )


class Protein(Document, Entity):
    id: UUID = UUIDField(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: ProteinName = ValueObjectStringField(required=True, factory=ProteinName)
    content = BinaryField(required=True)

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
            raise NoLabsException(ErrorCodes.invalid_protein_id)
        if not name:
            raise NoLabsException(ErrorCodes.invalid_protein_name)
        if not experiment:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)
        if not content:
            raise NoLabsException(ErrorCodes.invalid_protein_content)

        if isinstance(content, str):
            content = content.encode('utf-8')

        super().__init__(id=id.value if isinstance(id, ProteinId) else id,
                         experiment=experiment,
                         name=name,
                         content=content, *args, **kwargs)

        EventDispatcher.raise_event(ProteinCreatedEvent(self))

    @property
    def iid(self) -> ProteinId:
        return ProteinId(self.id)

    def set_name(self, name: ProteinName):
        if not name:
            raise NoLabsException(ErrorCodes.invalid_protein_name)
        self.name = name

    def set_content(self, content: Union[bytes, str]):
        if not content:
            raise NoLabsException(ErrorCodes.invalid_protein_content)

        if isinstance(content, str):
            content = content.encode('utf-8')

        self.content = content

    def get_content(self) -> str | None:
        if self.content:
            return self.content.decode('utf-8')

        return None

    @classmethod
    def create(cls, experiment: Experiment, name: ProteinName, content: Union[bytes, str]):
        proteins = cls.objects.filter(name=name.value, experiment=experiment)
        if proteins:
            protein: Protein = proteins[0]
            protein.set_content(content)
            protein.set_name(name)
            return protein

        return Protein(
            ProteinId(uuid.uuid4()),
            experiment=experiment,
            name=name,
            content=content
        )


class Job(Document, Entity):
    id: UUID = UUIDField(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: JobName = ValueObjectStringField(required=True, factory=JobName)
    created_at: datetime.datetime = DateTimeField(default=datetime.datetime.utcnow)

    meta = {
        'allow_inheritance': True
    }

    def __init__(self, id: JobId, name: JobName, experiment: Experiment, *args, **kwargs):
        if not id:
            raise NoLabsException(ErrorCodes.invalid_job_id)
        if not name:
            raise NoLabsException(ErrorCodes.invalid_job_name)
        if not experiment:
            raise NoLabsException(ErrorCodes.invalid_experiment_id)

        super().__init__(id=id.value if isinstance(id, JobId) else id, name=name, experiment=experiment, *args, **kwargs)

    def set_name(self, name: JobName):
        if not name:
            raise NoLabsException(ErrorCodes.invalid_job_name)

        self.name = name

    @property
    def iid(self) -> JobId:
        return JobId(self.id)


# region events
class AminoAcidCreatedEvent(DomainEvent):
    amino_acid: AminoAcid

    def __init__(self, amino_acid: AminoAcid):
        self.amino_acid = amino_acid


class ProteinCreatedEvent(DomainEvent):
    protein: Protein

    def __init__(self, protein: Protein):
        self.protein = protein

# endregion
