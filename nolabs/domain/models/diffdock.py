__all__ = ["DiffDockBindingJob"]

from typing import List, Optional

from mongoengine import CASCADE, PULL, IntField, ListField, ReferenceField, StringField

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Job, JobInputError, Ligand, Protein


class DiffDockBindingJob(Job):
    # region Inputs

    protein: Protein = ReferenceField(
        Protein, reverse_delete_rule=CASCADE, required=True
    )
    ligand: Ligand = ReferenceField(Ligand, reverse_delete_rule=CASCADE, required=True)
    samples_per_complex: int = IntField(default=2, required=False)

    # endregion

    celery_task_id: Optional[str] = StringField()

    complexes: List[Protein] = ListField(
        ReferenceField(Protein, required=True, reverse_delete_rule=PULL)
    )

    def set_task_id(self, task_id: str):
        self.celery_task_id = task_id

    def set_input(self, protein: Protein, ligand: Ligand, samples_per_complex: int = 2):
        self.complexes = []

        self.protein = protein
        self.ligand = ligand
        self.samples_per_complex = samples_per_complex

        self.input_errors(throw=True)

    def result_valid(self) -> bool:
        return not not self.complexes

    def input_valid(self) -> bool:
        return bool(self.protein and self.ligand)

    def set_result(self, complexes: List[Protein]):
        if not complexes:
            raise NoLabsException(ErrorCodes.invalid_job_result, "Result is empty")

        self.complexes = complexes

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
