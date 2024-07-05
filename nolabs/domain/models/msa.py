__all__ = [
    'MsaGenerationJob'
]

import datetime
from typing import List

from mongoengine import ReferenceField, CASCADE, BinaryField

from nolabs.domain.models.common import Job, Protein, JobInputError
from nolabs.exceptions import NoLabsException, ErrorCodes


class MsaGenerationJob(Job):
    # region Inputs

    protein: Protein | None = ReferenceField(Protein, required=True, reverse_delete_rule=CASCADE)
    msa: bytes | None = BinaryField(required=False)

    # endregion

    def set_input(self, protein: Protein):
        self.protein = protein

        self.input_errors(throw=True)

        self.inputs_updated_at = datetime.datetime.utcnow()

    def result_valid(self) -> bool:
        return not not self.msa

    def set_result(self,
                   protein: Protein,
                   msa: bytes | str):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty)

        if not msa:
            raise NoLabsException(ErrorCodes.invalid_msa)

        if isinstance(msa, str):
            msa = msa.encode('utf-8')

        self.msa = msa

    def _input_errors(self) -> List[JobInputError]:
        errors = []

        if not self.protein:
            errors.append([
                JobInputError(
                    message='Protein is not defined',
                    error_code=ErrorCodes.protein_is_undefined
                )
            ])

        if self.protein and not self.protein.fasta_content:
            errors.append([
                JobInputError(
                    message='Protein does not have fasta content',
                    error_code=ErrorCodes.protein_fasta_is_empty
                )
            ])

        return errors
