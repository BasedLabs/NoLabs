__all__ = [
    'LocalisationJob',
    'LocalisationProbability',
]

from typing import List, Callable, Awaitable

from mongoengine import EmbeddedDocument, ReferenceField, FloatField, ListField, EmbeddedDocumentListField, CASCADE, \
    PULL

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, JobType, AminoAcid, AminoAcidId
from nolabs.refined.infrastructure.mongo_fields import ValueObjectUUIDField
from nolabs.seedwork.domain.value_objects import ValueObject, ValueObjectUUID


class LocalisationProbability(EmbeddedDocument, ValueObject):
    amino_acid: AminoAcid = ReferenceField(AminoAcid, reverse_delete_rule=CASCADE)
    cytosolic_proteins: float = FloatField(required=True)
    mitochondrial_proteins: float = FloatField(required=True)
    nuclear_proteins: float = FloatField(required=True)
    other_proteins: float = FloatField(required=True)
    extracellular_secreted_proteins: float = FloatField(required=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        values = [
            self.cytosolic_proteins,
            self.mitochondrial_proteins,
            self.nuclear_proteins,
            self.other_proteins,
            self.extracellular_secreted_proteins
        ]

        for value in values:
            if not value:
                NoLabsException.throw(ErrorCodes.invalid_protein_location_probability)
            if value < 0 or value > 1.0:
                NoLabsException.throw(ErrorCodes.invalid_protein_location_probability)


class LocalisationJob(Job):
    amino_acids: List[AminoAcid] = ListField(ReferenceField(AminoAcid, required=False, reverse_delete_rule=PULL), required=False)
    result: List[LocalisationProbability] = EmbeddedDocumentListField(LocalisationProbability, required=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self.job_type != JobType.LOCALISATION:
            raise NoLabsException('Unsupported job type', ErrorCodes.invalid_job_type)

    def set_amino_acids(self, amino_acids: List[AminoAcid]):
        if not amino_acids:
            raise NoLabsException('Cannot set empty amino acids', ErrorCodes.invalid_job_input)

        self.amino_acids = amino_acids

    def set_result(self, result: List[LocalisationProbability]):
        if not self.amino_acids:
            raise NoLabsException('Cannot set result on empty inputs', ErrorCodes.invalid_job_input)

        if not result:
            raise NoLabsException('Cannot set empty result', ErrorCodes.invalid_job_result)

        self.result = result
