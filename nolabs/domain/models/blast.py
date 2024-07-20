__all__ = [
    'DiffDockBindingJob'
]

import datetime
from typing import List
from uuid import UUID

from mongoengine import ReferenceField, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, CASCADE, IntField, BinaryField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import Job, Protein, JobInputError


class BlastJobResult(EmbeddedDocument):
    complex_id: UUID = UUIDField(required=True)
    sdf_content: bytes = BinaryField(required=True)
    minimized_affinity: float = FloatField(required=False)
    scored_affinity: float = FloatField(required=False)
    confidence: float = FloatField(required=False)

    @staticmethod
    def create(complex_id: UUID,
               sdf_content: bytes | str,
               minimized_affinity: float,
               scored_affinity: float,
               confidence: float
               ):
        if isinstance(sdf_content, str):
            sdf_content = sdf_content.encode('utf-8')

        return BlastJobResult(
            complex_id=complex_id,
            sdf_content=sdf_content,
            minimized_affinity=minimized_affinity,
            scored_affinity=scored_affinity,
            confidence=confidence
        )


class BlastJob(Job):
    # region Inputs

    protein: Protein = ReferenceField(Protein, reverse_delete_rule=CASCADE, required=True)

    # endregion

    result: List[BlastJobResult] = EmbeddedDocumentListField(BlastJobResult)

    def set_input(self,
                  protein: Protein,
                  samples_per_complex: int = 1
                  ):
        self.result = []

        self.protein = protein
        self.samples_per_complex = samples_per_complex

        self.input_errors(throw=True)

        self.inputs_updated_at = datetime.datetime.utcnow()

    def result_valid(self) -> bool:
        return not not self.result

    def input_valid(self) -> bool:
        return bool(self.protein and self.ligand)

    def set_result(self,
                   protein: Protein,
                   result: List[BlastJobResult]):
        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Cannot set a result on a job without inputs')

        if not protein:
            raise NoLabsException(ErrorCodes.protein_is_undefined, 'Provided protein is undefined')

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs, 'Protein does not match with job inputs')

        if not result:
            raise NoLabsException(ErrorCodes.invalid_job_result, 'Result is empty')

        self.result = result

    def _input_errors(self) -> List[JobInputError]:
        errors = []

        if not self.protein:
            errors.append(
                JobInputError(
                    message='Protein is undefined',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        if self.protein and not self.protein.pdb_content:
            errors.append(
                JobInputError(
                    message='Protein pdb content is undefined',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        return errors
