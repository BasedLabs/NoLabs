__all__ = ["FoldingJobResult", "FoldingJob"]

from datetime import datetime
from enum import Enum
from typing import List
from uuid import UUID

from domain.exceptions import ErrorCodes, NoLabsException
from mongoengine import (CASCADE, BinaryField, EmbeddedDocument,
                         EmbeddedDocumentField, EnumField, ReferenceField,
                         UUIDField)

from nolabs.domain.models.common import Job, JobInputError, Protein


class FoldingBackendEnum(str, Enum):
    esmfold = "esmfold"
    esmfold_light = "esmfold_light"
    rosettafold = "rosettafold"


class FoldingJobResult(EmbeddedDocument):
    """
    Not a domain object
    Currently used just to keep job outputs
    """

    protein_id: UUID = UUIDField(required=True)
    pdb_content: bytes = BinaryField(required=True)


class FoldingJob(Job):
    # region Inputs

    protein: Protein = ReferenceField(
        Protein, required=True, reverse_delete_rule=CASCADE
    )
    backend: FoldingBackendEnum = EnumField(FoldingBackendEnum, required=True)

    # endregion

    folding: FoldingJobResult = EmbeddedDocumentField(FoldingJobResult)

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

        self.inputs_updated_at = datetime.utcnow()

    def result_valid(self) -> bool:
        return not not self.folding

    def set_result(self, protein: Protein, pdb: str | bytes):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if not pdb:
            raise NoLabsException(ErrorCodes.invalid_job_result)

        self.folding = None

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        self.folding = FoldingJobResult(
            protein_id=protein.iid.value,
            pdb_content=pdb if isinstance(pdb, bytes) else pdb.encode("utf-8"),
        )

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
