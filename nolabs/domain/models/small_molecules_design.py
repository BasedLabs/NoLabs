__all__ = ["SmallMoleculesDesignJob"]

from typing import List, Optional

from mongoengine import (
    CASCADE,
    PULL,
    FloatField,
    IntField,
    ListField,
    ReferenceField,
    StringField,
)

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Job, JobInputError, Ligand, Protein


class SmallMoleculesDesignJob(Job):
    # region Inputs

    protein: Protein = ReferenceField(
        Protein, required=True, reverse_delete_rule=CASCADE
    )

    # endregion

    center_x: float = FloatField(default=0.0)
    center_y: float = FloatField(default=0.0)
    center_z: float = FloatField(default=0.0)
    size_x: float = FloatField(default=5.0)
    size_y: float = FloatField(default=5.0)
    size_z: float = FloatField(default=5.0)
    batch_size: int = IntField(default=128)
    minscore: float = FloatField(default=0.4)
    epochs: int = IntField(default=50)

    ligands: List[Ligand] = ListField(
        ReferenceField(Ligand, required=False, reverse_delete_rule=PULL)
    )

    celery_task_id: Optional[str] = StringField()

    def result_valid(self) -> bool:
        return not not self.ligands

    def set_protein(self, protein: Protein):
        if not protein:
            raise NoLabsException(ErrorCodes.protein_is_undefined)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_not_found)

        self.protein = protein

    def set_inputs(
        self,
        protein: Protein,
        center_x: float,
        center_y: float,
        center_z: float,
        size_x: float,
        size_y: float,
        size_z: float,
        batch_size: int,
        minscore: float,
        epochs: int,
        throw: bool = True,
    ):
        if (
            self.center_x != center_x
            or self.center_y != center_y
            or self.center_z != center_z
            or self.size_x != size_x
            or self.size_y != size_y
            or self.size_z != size_z
            or self.batch_size != batch_size
            or self.minscore != minscore
            or self.epochs != epochs
            or not self.protein
            or not self.protein.id != protein.id
        ):
            self.ligands = []
            self.protein = protein
            self.center_x = center_x
            self.center_y = center_y
            self.center_z = center_z
            self.size_x = size_x
            self.size_y = size_y
            self.size_z = size_z
            self.batch_size = batch_size
            self.minscore = minscore
            self.epochs = epochs

            self.processing_required = True

        if throw:
            self.input_errors(throw=True)

    def set_result(self, ligands: List[Ligand]):
        if not ligands:
            raise NoLabsException(ErrorCodes.small_molecules_design_empty_output)

        self.ligands = ligands
        self.processing_required = False

    def set_task_id(self, task_id):
        self.celery_task_id = task_id

    def _input_errors(self) -> List[JobInputError]:
        errors = []

        if not self.protein:
            errors.append(
                JobInputError(
                    message="Protein is undefined",
                    error_code=ErrorCodes.protein_is_undefined,
                )
            )

        if not self.protein.pdb_content:
            errors.append(
                JobInputError(
                    message="Protein pdb content is empty",
                    error_code=ErrorCodes.protein_is_undefined,
                )
            )

        if not self.center_x:
            errors.append(
                JobInputError(
                    message="Center x of binding box is undefined or 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.center_y:
            errors.append(
                JobInputError(
                    message="Center н of binding box is undefined or 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.center_z:
            errors.append(
                JobInputError(
                    message="Center z of binding box is undefined or 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.size_x:
            errors.append(
                JobInputError(
                    message="Size x of binding box is undefined or 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.size_y:
            errors.append(
                JobInputError(
                    message="Size y of binding box is undefined or 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.size_y:
            errors.append(
                JobInputError(
                    message="Size y of binding box is undefined or 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.batch_size:
            errors.append(
                JobInputError(
                    message="Batch size is undefined or 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.minscore:
            errors.append(
                JobInputError(
                    message="Min molecule acceptance score is undefined or 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.epochs or self.epochs < 0:
            errors.append(
                JobInputError(
                    message="Number of epochs are less or equal to 0",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        return errors
