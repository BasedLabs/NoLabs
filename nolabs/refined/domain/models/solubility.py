__all__ = [
    'SolubilityJobResult',
    'SolubilityJob'
]

from enum import Enum
from typing import List
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, BinaryField, StringField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, LocalisationProbability, SolubleProbability


class SolubilityJobResult(EmbeddedDocument):
    protein_id: UUID = UUIDField(required=True)
    soluble_probability: float = FloatField(required=True)


class SolubilityJob(Job):
    proteins: List[Protein] = ListField(ReferenceField(Protein, required=False, reverse_delete_rule=PULL))
    results: List[SolubilityJobResult] = EmbeddedDocumentListField(SolubilityJobResult)

    def set_proteins(self, proteins: List[Protein]):
        if not proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        self.proteins = proteins

    def clear_result(self):
        self.results = []

    def set_result(self, protein: Protein, soluble_probability: SolubleProbability):
        if not soluble_probability:
            raise NoLabsException(ErrorCodes.invalid_solubility_probability)

        if not self.proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein not in self.proteins:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.fasta_content:
            raise NoLabsException(ErrorCodes.protein_amino_acid_sequence_not_found)

        existing_result = [res for res in self.results if res.protein_id == protein.id]
        if existing_result:
            soluble_probability_existing = existing_result[0]
            soluble_probability_existing.soluble_probability = soluble_probability.value
        else:
            result = SolubilityJobResult(
                protein_id=protein.id,
                soluble_probability=soluble_probability.value
            )
            self.results.append(result)
