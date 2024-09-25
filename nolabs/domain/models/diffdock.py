__all__ = ["DiffDockBindingJob"]

import datetime
import uuid
from typing import List
from uuid import UUID

from mongoengine import (CASCADE, BinaryField, EmbeddedDocument,
                         EmbeddedDocumentListField, FloatField, IntField,
                         ReferenceField, UUIDField)

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Job, JobInputError, Ligand, Protein


class DiffDockJobResult(EmbeddedDocument):
    complex_id: UUID = UUIDField(required=True)
    sdf_content: bytes = BinaryField(required=True)
    minimized_affinity: float = FloatField(required=False)
    scored_affinity: float = FloatField(required=False)
    confidence: float = FloatField(required=False)

    @staticmethod
    def create(
        complex_id: UUID,
        sdf_content: bytes | str,
        minimized_affinity: float,
        scored_affinity: float,
        confidence: float,
    ):
        if isinstance(sdf_content, str):
            sdf_content = sdf_content.encode("utf-8")

        return DiffDockJobResult(
            complex_id=complex_id,
            sdf_content=sdf_content,
            minimized_affinity=minimized_affinity,
            scored_affinity=scored_affinity,
            confidence=confidence,
        )


class DiffDockBindingJob(Job):
    # region Inputs

    protein: Protein = ReferenceField(
        Protein, reverse_delete_rule=CASCADE, required=True
    )
    ligand: Ligand = ReferenceField(Ligand, reverse_delete_rule=CASCADE, required=True)
    samples_per_complex: int = IntField(default=2, required=False)

    # endregion

    celery_task_id: uuid.UUID = UUIDField()

    result: List[DiffDockJobResult] = EmbeddedDocumentListField(DiffDockJobResult)

    def set_input(self, protein: Protein, ligand: Ligand, samples_per_complex: int = 1):
        self.result = []

        self.protein = protein
        self.ligand = ligand
        self.samples_per_complex = samples_per_complex

        self.input_errors(throw=True)

        self.inputs_updated_at = datetime.datetime.utcnow()

    def result_valid(self) -> bool:
        return not not self.result

    def input_valid(self) -> bool:
        return bool(self.protein and self.ligand)

    def set_result(self, result: List[DiffDockJobResult]):
        if not result:
            raise NoLabsException(ErrorCodes.invalid_job_result, "Result is empty")

        self.result = result

    def _input_errors(self) -> List[JobInputError]:
        errors = []

        if not self.protein:
            errors.append(
                JobInputError(
                    message="Protein is undefined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if self.protein and not self.protein.pdb_content:
            errors.append(
                JobInputError(
                    message="Protein pdb content is undefined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.ligand:
            errors.append(
                JobInputError(
                    message="Ligand is undefined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if self.ligand and not self.ligand.sdf_content:
            errors.append(
                JobInputError(
                    message="Ligand sdf content is undefined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if self.samples_per_complex <= 0:
            errors.append(
                JobInputError(
                    message="Samples per complex is undefined or <= 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        return errors
