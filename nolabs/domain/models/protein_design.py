__all__ = ["ProteinDesignJob"]

import datetime
from typing import List

from mongoengine import CASCADE, PULL, IntField, ListField, ReferenceField, StringField

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Job, JobInputError, Protein


class ProteinDesignJob(Job):
    # region Inputs

    protein: Protein | None = ReferenceField(
        Protein, required=True, reverse_delete_rule=CASCADE
    )
    binders: List[Protein] = ListField(
        ReferenceField(Protein, required=True, reverse_delete_rule=PULL)
    )

    # endregion

    contig: str = StringField(required=False)
    number_of_designs: int = IntField(required=False, default=2)
    hotspots: str = StringField(required=False)
    timesteps: int = IntField(required=False, default=50)

    def set_input(
        self,
        protein: Protein,
        contig: str,
        number_of_designs: int,
        hotspots: str,
        timesteps: int,
    ):
        self.binders = []
        self.contig = contig
        self.number_of_designs = number_of_designs
        self.hotspots = hotspots
        self.timesteps = timesteps
        self.protein = protein

        self.input_errors(throw=True)

        self.inputs_updated_at = datetime.datetime.utcnow()

    def set_protein(self, protein: Protein):
        if not protein:
            raise NoLabsException(ErrorCodes.protein_is_undefined)

        if not protein.pdb_content:
            raise NoLabsException

        self.protein = protein
        self.inputs_updated_at = datetime.datetime.utcnow()

    def result_valid(self) -> bool:
        return not not self.binders

    def set_result(self, protein: Protein, binders: List[Protein]):
        self.binders = []

        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty)

        self.binders = binders

    def _input_errors(self) -> List[JobInputError]:
        errors = []

        if not self.protein:
            errors.append(
                JobInputError(
                    message="Protein is not defined",
                    error_code=ErrorCodes.protein_is_undefined,
                )
            )

        if self.protein and not self.protein.pdb_content:
            errors.append(
                JobInputError(
                    message="Protein does not have pdb content",
                    error_code=ErrorCodes.protein_pdb_is_empty,
                )
            )

        if not self.contig:
            errors.append(
                JobInputError(
                    message="Contig is undefined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.number_of_designs or self.number_of_designs <= 0:
            errors.append(
                JobInputError(
                    message="Number of designs cannot be <= 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.timesteps or self.timesteps <= 0:
            errors.append(
                JobInputError(
                    message="Timesteps cannot be <= 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        return errors
