__all__ = [
    'PocketPredictionJob'
]

from typing import List

from mongoengine import ReferenceField, ListField, CASCADE, IntField

from nolabs.domain.models.common import Job, Protein, JobInputError
from nolabs.exceptions import NoLabsException, ErrorCodes


class PocketPredictionJob(Job):
    # region Inputs

    protein: Protein | None = ReferenceField(Protein, required=True, reverse_delete_rule=CASCADE)

    # endregion

    pocket_ids: List[int] = ListField(IntField())

    def set_input(self,
                  protein: Protein
                  ):
        self.pocket_ids = []
        self.protein = protein

        self.input_errors(throw=True)

    def result_valid(self) -> bool:
        return not not self.pocket_ids

    def set_result(self,
                   protein: Protein,
                   pocket_ids: List[int]):
        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty)

        self.pocket_ids = pocket_ids

    def _input_errors(self) -> List[JobInputError]:
        errors = []

        if not self.protein:
            errors.append(JobInputError(message='Protein is not defined',
                                        error_code=ErrorCodes.protein_is_undefined))

        if not self.protein.pdb_content:
            errors.append(JobInputError(message='Protein pdb content is empty',
                                  error_code=ErrorCodes.protein_pdb_is_empty))

        return errors
