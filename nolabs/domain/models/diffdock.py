__all__ = [
    'DiffDockBindingJob'
]

import datetime
from typing import List
from uuid import UUID

from mongoengine import ReferenceField, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, CASCADE, IntField, BinaryField

from domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import Job, Protein, Ligand, JobInputError


class DiffDockJobResult(EmbeddedDocument):
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

        return DiffDockJobResult(
            complex_id=complex_id,
            sdf_content=sdf_content,
            minimized_affinity=minimized_affinity,
            scored_affinity=scored_affinity,
            confidence=confidence
        )


class DiffDockBindingJob(Job):
    # region Inputs

    protein: Protein = ReferenceField(Protein, reverse_delete_rule=CASCADE, required=True)
    ligand: Ligand = ReferenceField(Ligand, reverse_delete_rule=CASCADE, required=True)
    samples_per_complex: int = IntField(default=2, required=False)

    # endregion

    result: List[DiffDockJobResult] = EmbeddedDocumentListField(DiffDockJobResult)

    def set_input(self,
                  protein: Protein,
                  ligand: Ligand,
                  samples_per_complex: int = 1
                  ):
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

    def set_result(self,
                   protein: Protein,
                   ligand: Ligand,
                   result: List[DiffDockJobResult]):
        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Cannot set a result on a job without inputs')

        if not protein:
            raise NoLabsException(ErrorCodes.protein_is_undefined, 'Provided protein is undefined')

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs, 'Protein does not match with job inputs')

        if not ligand:
            raise NoLabsException(ErrorCodes.ligand_is_undefined, 'Provided ligand is undefined')

        if ligand != self.ligand:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Some or all of ligands were not a part of job input')

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

        if not self.ligand:
            errors.append(
                JobInputError(
                    message='Ligand is undefined',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        if self.ligand and not self.ligand.sdf_content:
            errors.append(
                JobInputError(
                    message='Ligand sdf content is undefined',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        if self.samples_per_complex <= 0:
            errors.append(
                JobInputError(
                    message='Samples per complex is undefined or <= 0',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        return errors
