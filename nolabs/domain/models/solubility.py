__all__ = ["SolubilityJobResult", "SolubilityJob"]

from datetime import datetime
from typing import List, Tuple
from uuid import UUID

from domain.exceptions import ErrorCodes, NoLabsException
from mongoengine import (PULL, EmbeddedDocument, EmbeddedDocumentListField,
                         FloatField, ListField, ReferenceField, UUIDField)

from nolabs.domain.models.common import Job, JobInputError, Protein


class SolubilityJobResult(EmbeddedDocument):
    protein_id: UUID = UUIDField(required=True)
    soluble_probability: float = FloatField(required=True)


class SolubilityJob(Job):
    # region Inputs

    proteins: List[Protein] = ListField(
        ReferenceField(Protein, required=True, reverse_delete_rule=PULL)
    )

    # endregion

    results: List[SolubilityJobResult] = EmbeddedDocumentListField(SolubilityJobResult)

    def set_input(self, proteins: List[Protein]):
        self.proteins = proteins

        self.input_errors(throw=True)

        self.inputs_updated_at = datetime.utcnow()

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

            self.results.append(
                SolubilityJobResult(
                    protein_id=protein.iid.value, soluble_probability=prob
                )
            )

    def _input_errors(self) -> List[JobInputError]:
        if not self.proteins:
            return [
                JobInputError(
                    message="No proteins provided",
                    error_code=ErrorCodes.protein_is_undefined,
                )
            ]

        for protein in self.proteins:
            if not protein.fasta_content:
                return [
                    JobInputError(
                        message="Protein fasta is undefined",
                        error_code=ErrorCodes.protein_fasta_is_empty,
                    )
                ]

        return []
