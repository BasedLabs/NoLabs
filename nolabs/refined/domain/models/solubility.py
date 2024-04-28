from typing import List

from mongoengine import FloatField, EmbeddedDocument, ReferenceField, ListField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models import AminoAcid, Job
from nolabs.seedwork.domain.value_objects import ValueObject


class SolubleProbability(EmbeddedDocument, ValueObject):
    amino_acid: AminoAcid = ReferenceField(AminoAcid)
    value: float = FloatField(required=True)

    def clean(self):
        if self.value < 0 or self.value > 1.0:
            NoLabsException.throw(ErrorCodes.invalid_protein_location_probability)


class SolubilityJob(Job):
    amino_acids: List[AminoAcid] = ListField(ReferenceField(AminoAcid), required=False)
    result: List[SolubleProbability] = ListField(EmbeddedDocument(SolubleProbability), required=False)

    def clean(self):
        super().clean()

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

    async def start(self, predictor: Callable[[List[AminoAcid]], Awaitable[List[LocalisationProbability]]]):
        if not self.amino_acids:
            raise NoLabsException('Cannot start job with empty inputs', ErrorCodes.invalid_job_input)

        return predictor(self.amino_acids)