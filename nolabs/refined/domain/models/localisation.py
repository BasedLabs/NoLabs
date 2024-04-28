__all__ = [
    'LocalisationJob',
    'LocalisationProbability',
]

import datetime
from typing import List

from mongoengine import EmbeddedDocument, ReferenceField, FloatField, ListField, EmbeddedDocumentListField, CASCADE, \
    PULL

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, JobType, AminoAcid, JobId, JobName
from nolabs.refined.domain.models.experiment import Experiment
from nolabs.seedwork.domain.value_objects import ValueObject


class LocalisationProbability(EmbeddedDocument, ValueObject):
    amino_acid: AminoAcid = ReferenceField(AminoAcid, reverse_delete_rule=CASCADE)
    cytosolic_proteins: float = FloatField(required=True)
    mitochondrial_proteins: float = FloatField(required=True)
    nuclear_proteins: float = FloatField(required=True)
    other_proteins: float = FloatField(required=True)
    extracellular_secreted_proteins: float = FloatField(required=True)

    def __init__(self, amino_acid: AminoAcid, cytosolic_proteins: float, mitochondrial_proteins: float,
                 nuclear_proteins: float, other_proteins: float, extracellular_secreted_proteins: float, *args,
                 **kwargs):
        values = [
            cytosolic_proteins,
            mitochondrial_proteins,
            nuclear_proteins,
            other_proteins,
            extracellular_secreted_proteins
        ]

        for value in values:
            if not value:
                NoLabsException.throw(ErrorCodes.invalid_protein_location_probability)
            if value < 0 or value > 1.0:
                NoLabsException.throw(ErrorCodes.invalid_protein_location_probability)

        if not amino_acid:
            raise NoLabsException.throw(ErrorCodes.invalid_aa_id)

        super().__init__(
            amino_acid=amino_acid,
            cytosolic_proteins=cytosolic_proteins,
            mitochondrial_proteins=mitochondrial_proteins,
            nuclear_proteins=nuclear_proteins,
            other_proteins=other_proteins,
            extracellular_secreted_proteins=extracellular_secreted_proteins,
            *args, **kwargs)

    def clean(self):
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
    amino_acids: List[AminoAcid] = ListField(ReferenceField(AminoAcid, required=False, reverse_delete_rule=PULL),
                                             required=False)
    result: List[LocalisationProbability] = EmbeddedDocumentListField(LocalisationProbability, required=False)

    def __init__(self, id: JobId, name: JobName, experiment: Experiment, *args, **values):
        if not id:
            raise NoLabsException('Job id is invalid', ErrorCodes.invalid_job_id)
        if not name:
            raise NoLabsException('Job name is invalid', ErrorCodes.invalid_job_name)
        if not experiment:
            raise NoLabsException('Job must be in experiment', ErrorCodes.invalid_experiment_id)

        super().__init__(id=id, name=name, experiment=experiment, *args, **values)

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
