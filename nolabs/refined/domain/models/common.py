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
from enum import Enum
from functools import cached_property
from typing import Union
from uuid import UUID

from mongoengine import EnumField, DateTimeField, FileField, Document, ReferenceField, CASCADE, EmbeddedDocument, \
    FloatField, EmbeddedDocumentField, BinaryField, fields
from pydantic.dataclasses import dataclass

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.event_dispatcher import EventDispatcher
from nolabs.refined.infrastructure.mongo_fields import ValueObjectUUIDField, ValueObjectStringField
from nolabs.seedwork.domain.entities import Entity
from nolabs.seedwork.domain.events import DomainEvent
from nolabs.seedwork.domain.value_objects import ValueObject, ValueObjectString, ValueObjectUUID


@dataclass
class ExperimentId(ValueObjectUUID):
    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_experiment_id)

    def __eq__(self, other):
        if not isinstance(other, ProteinId):
            return False

        return self.value == other.value

    class Field(fields.BaseField):
        def to_mongo(self, value):
            if isinstance(value, ExperimentId):
                return value.value
            else:
                return value

        def to_python(self, value):
            if isinstance(value, UUID):
                return ExperimentId(value)
            return value


@dataclass
class ExperimentName(ValueObjectString):
    def __post_init__(self):
        if not self.value:
            raise NoLabsException('Experiment name cannot be empty', ErrorCodes.invalid_experiment_name)

        if len(self.value) < 1 or len(self.value) > 1000:
            raise NoLabsException('Length of experiment name must be between 1 and 1000',
                                  ErrorCodes.invalid_experiment_name)

    class Field(fields.BaseField):
        def to_mongo(self, value):
            if isinstance(value, ExperimentName):
                return value.value
            else:
                return value

        def to_python(self, value):
            if isinstance(value, str):
                return ExperimentName(value)
            return value


class Experiment(Document, Entity):
    id: ExperimentId = ExperimentId.Field(db_field='_id', primary_key=True, required=True)
    name: ExperimentName = ExperimentName.Field(required=True)
    created_at: datetime.datetime = DateTimeField(default=datetime.datetime.utcnow)

    def __init__(self, id: ExperimentId, name: ExperimentName, created_at: datetime.datetime | None = None, *args, **kwargs):
        if not id:
            raise NoLabsException('Experiment id is empty', ErrorCodes.invalid_experiment_id)

        if not name:
            raise NoLabsException('Experiment name is empty', ErrorCodes.invalid_experiment_name)

        created_at = created_at if created_at else datetime.datetime.now(tz=datetime.timezone.utc)

        super().__init__(id=id, name=name, created_at=created_at, *args, **kwargs)

    def set_name(self, name: ExperimentName):
        if not name:
            raise NoLabsException('Name cannot be empty', ErrorCodes.invalid_experiment_name)

        self.name = name

    def delete(self, signal_kwargs=None, **write_concern):
        self.__class__.objects.filter(id=self.id.value).delete()


@dataclass
class ProteinName(ValueObjectString):
    def __post_init__(self):
        if not self.value:
            raise NoLabsException('Protein name cannot be empty', ErrorCodes.invalid_protein_name)

        if len(self.value) < 1 or len(self.value) > 1000:
            raise NoLabsException('Length of protein name must be between 1 and 1000', ErrorCodes.invalid_protein_name)

    class Field(fields.BaseField):
        def to_mongo(self, value):
            if isinstance(value, ProteinName):
                return value.value
            else:
                return value

        def to_python(self, value):
            if isinstance(value, str):
                return ProteinName(value)
            return value


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

    class Field(fields.BaseField):
        def to_mongo(self, value):
            if isinstance(value, AminoAcidName):
                return value.value
            else:
                return value

        def to_python(self, value):
            if isinstance(value, str):
                return AminoAcidName(value)
            return value


@dataclass
class ProteinId(ValueObjectUUID):
    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_protein_id)

    def __eq__(self, other):
        if not isinstance(other, ProteinId):
            return False

        return self.value == other.value

    class Field(fields.BaseField):
        def to_mongo(self, value):
            if isinstance(value, ProteinId):
                return value.value
            else:
                return value

        def to_python(self, value):
            if isinstance(value, UUID):
                return ProteinId(value)
            return value


@dataclass
class AminoAcidId(ValueObjectUUID):
    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_aa_id)

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if not isinstance(other, AminoAcidId):
            return False

        return self.value == other.value

    class Field(fields.BaseField):
        def to_mongo(self, value):
            if isinstance(value, AminoAcidId):
                return value.value
            else:
                return value

        def to_python(self, value):
            if isinstance(value, UUID):
                return AminoAcidId(value)
            return value


@dataclass
class JobId(ValueObjectUUID):
    value: UUID

    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_job_id)

    def __eq__(self, other):
        if not isinstance(other, ProteinId):
            return False

        return self.value == other.value

    class Field(fields.BaseField):
        def to_mongo(self, value):
            if isinstance(value, JobId):
                return value.value
            else:
                return value

        def to_python(self, value):
            if isinstance(value, UUID):
                return JobId(value)
            return value


class JobName(ValueObjectString):
    def __post_init__(self):
        if not self.value:
            NoLabsException.throw(ErrorCodes.invalid_job_name)

        if len(self.value) < 1 or len(self.value) > 100:
            NoLabsException.throw(ErrorCodes.invalid_job_name)

    class Field(fields.BaseField):
        def to_mongo(self, value):
            if isinstance(value, JobName):
                return value.value
            else:
                return value

        def to_python(self, value):
            if isinstance(value, str):
                return JobName(value)
            return value


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
                NoLabsException.throw(ErrorCodes.invalid_protein_location_probability)
            if value < 0 or value > 1.0:
                NoLabsException.throw(ErrorCodes.invalid_protein_location_probability)

        super().__init__(
            cytosolic=cytosolic,
            mitochondrial=mitochondrial,
            nuclear=nuclear,
            other=other,
            extracellular=extracellular,
            *args,
            **kwargs)


class AminoAcid(Document, Entity):
    id: AminoAcidId = AminoAcidId.Field(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: AminoAcidName = AminoAcidName.Field(required=True)
    content = BinaryField(required=True)

    localisation: LocalisationProbability = EmbeddedDocumentField(LocalisationProbability, required=False)

    def __init__(self, id: AminoAcidId, experiment: Experiment, name: AminoAcidName, content: Union[bytes, str], *args,
                 **kwargs):
        if not id:
            raise NoLabsException('Amino acid id is invalid', ErrorCodes.invalid_aa_id)
        if not name:
            raise NoLabsException('Amino acid name is invalid', ErrorCodes.invalid_aa_name)
        if not content:
            raise NoLabsException('Amino acid content is invalid', ErrorCodes.invalid_aa_content)
        if not experiment:
            raise NoLabsException('Amino acid must be in experiment', ErrorCodes.invalid_experiment_id)

        if isinstance(content, str):
            content = content.encode('utf-8')

        super().__init__(id=id, experiment=experiment, name=name, content=content, *args, **kwargs)

        EventDispatcher.raise_event(AminoAcidCreatedEvent(self))

    def set_name(self, name: AminoAcidName):
        if not name:
            NoLabsException.throw(ErrorCodes.invalid_aa_name)
        self.name = name

    def set_content(self, content: Union[bytes, str]):
        if not content:
            raise NoLabsException('Amino acid content is invalid', ErrorCodes.invalid_aa_content)

        if isinstance(content, str):
            content = content.encode('utf-8')

        self.content = content

    def get_content(self) -> str | None:
        if self.content:
            return self.content.decode('utf-8')

        return None

    def set_localisation_probability(self, localisation: LocalisationProbability):
        if not localisation:
            raise NoLabsException('Localisation probability is empty', ErrorCodes.invalid_localisation_probability)

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

    def delete(self, signal_kwargs=None, **write_concern):
        self.__class__.objects.filter(id=self.id.value).delete()


class Protein(Document, Entity):
    id: ProteinId = ProteinId.Field(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: ProteinName = ProteinName.Field(required=True)
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
            raise NoLabsException('Protein id is invalid', ErrorCodes.invalid_protein_id)
        if not name:
            raise NoLabsException('Protein name is invalid', ErrorCodes.invalid_protein_name)
        if not experiment:
            raise NoLabsException('Protein must be in experiment', ErrorCodes.invalid_experiment_id)
        if not content:
            raise NoLabsException('Protein content is empty', ErrorCodes.invalid_protein_content)

        if isinstance(content, str):
            content = content.encode('utf-8')

        super().__init__(id=id, experiment=experiment, name=name, content=content, *args, **kwargs)

        EventDispatcher.raise_event(ProteinCreatedEvent(self))

    def set_name(self, name: ProteinName):
        if not name:
            NoLabsException.throw(ErrorCodes.invalid_protein_name)
        self.name = name

    def set_content(self, content: Union[bytes, str]):
        if not content:
            raise NoLabsException('Protein content is invalid', ErrorCodes.invalid_protein_content)

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

    def delete(self, signal_kwargs=None, **write_concern):
        self.__class__.objects.filter(id=self.id.value).delete()


class Job(Document, Entity):
    id: JobId = JobId.Field(db_field='_id', primary_key=True, required=True)
    experiment: Experiment = ReferenceField(Experiment, required=True, reverse_delete_rule=CASCADE)
    name: JobName = JobName.Field(required=True)
    created_at: datetime.datetime = DateTimeField(default=datetime.datetime.utcnow)

    meta = {
        'allow_inheritance': True
    }

    def delete(self, signal_kwargs=None, **write_concern):
        self.__class__.objects.filter(id=self.id.value).delete()


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
