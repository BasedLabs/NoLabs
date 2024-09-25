__all__ = ["FoldingJob"]

from enum import Enum
from typing import List

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from mongoengine import (CASCADE, BinaryField, EmbeddedDocument,
                         EmbeddedDocumentField, EnumField, ReferenceField,
                         UUIDField)

from nolabs.domain.models.common import Job, JobInputError, Protein


class FoldingBackendEnum(str, Enum):
    esmfold = "esmfold"
    esmfold_light = "esmfold_light"
    rosettafold = "rosettafold"


class FoldingJob(Job):
    protein: Protein = ReferenceField(
        Protein, required=True, reverse_delete_rule=CASCADE
    )
    folded_protein: Protein = ReferenceField(Protein, required=False)
    backend: FoldingBackendEnum = EnumField(FoldingBackendEnum, required=True)

    def set_inputs(self, protein: Protein, backend: FoldingBackendEnum):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input, "Protein was not found")

        if not Protein.objects.with_id(protein.iid.value):
            raise NoLabsException(
                ErrorCodes.invalid_job_input,
                "Protein was not found in list of proteins",
            )
        if not protein.fasta_content:
            raise NoLabsException(
                ErrorCodes.protein_fasta_is_empty,
                "Cannot run folding on a protein without fasta content",
            )

        self.backend = backend
        self.protein = protein
        self.processing_required = True

    def result_valid(self) -> bool:
        return not not self.folded_protein

    def set_result(self, protein: Protein):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        self.folded_protein = protein
        self.processing_required = False

    def _input_errors(self) -> List[JobInputError]:
        if not self.protein:
            return [
                JobInputError(
                    message="Protein is undefined",
                    error_code=ErrorCodes.protein_is_undefined,
                )
            ]

        if not self.protein.fasta_content:
            return [
                JobInputError(
                    message="Protein fasta content is undefined",
                    error_code=ErrorCodes.protein_fasta_is_empty,
                )
            ]

        return []
