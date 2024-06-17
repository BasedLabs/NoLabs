__all__ = [
    'SolubilityJobResult',
    'SolubilityJob'
]

from enum import Enum
from typing import List, Tuple
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, BinaryField, StringField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, LocalisationProbability, SolubleProbability


class SolubilityJobResult(EmbeddedDocument):
    protein_id: UUID = UUIDField(required=True)
    soluble_probability: float = FloatField(required=True)


class SolubilityJob(Job):
    # region Inputs

    proteins: List[Protein] = ListField(ReferenceField(Protein, required=True, reverse_delete_rule=PULL))

    # endregion

    results: List[SolubilityJobResult] = EmbeddedDocumentListField(SolubilityJobResult)

    def set_proteins(self, proteins: List[Protein]):
        if not proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        for protein in proteins:
            if not protein.fasta_content:
                raise NoLabsException(ErrorCodes.protein_fasta_is_empty)

        self.proteins = proteins

    def result_valid(self) -> bool:
        return not not self.results

    def set_result(self, result: List[Tuple[Protein, float]]):
        if not self.proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        for protein in self.proteins:
            if not protein.fasta_content:
                raise NoLabsException(ErrorCodes.protein_fasta_is_empty)

        if not result:
            raise NoLabsException(ErrorCodes.invalid_job_result)

        self.results = []

        for protein, prob in result:
            if not [p for p in self.proteins if p.iid == protein.iid]:
                raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

            self.results.append(SolubilityJobResult(
                protein_id=protein.iid.value,
                soluble_probability=prob
            ))
