__all__ = [
    'FoldingJobResult',
    'FoldingJob'
]

from enum import Enum
from typing import List, Tuple
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, EmbeddedDocumentListField, \
    UUIDField, BinaryField, EnumField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import Job, Protein


class FoldingBackendEnum(str, Enum):
    esmfold = 'esmfold'
    esmfold_light = 'esmfold_light'
    rosettafold = 'rosettafold'


class FoldingJobResult(EmbeddedDocument):
    """
    Not a domain object
    Currently used just to keep job outputs
    """
    protein_id: UUID = UUIDField(required=True)
    pdb_content: bytes = BinaryField(required=True)


class FoldingJob(Job):
    # region Inputs

    proteins: List[Protein] = ListField(ReferenceField(Protein, required=True, reverse_delete_rule=PULL))
    backend: FoldingBackendEnum = EnumField(FoldingBackendEnum, required=True)

    # endregion

    foldings: List[FoldingJobResult] = EmbeddedDocumentListField(FoldingJobResult)

    def set_inputs(self, proteins: List[Protein], backend: FoldingBackendEnum):
        self.foldings = []

        if not proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Protein was not found')

        for protein in proteins:
            if not Protein.objects.with_id(protein.iid.value):
                raise NoLabsException(ErrorCodes.invalid_job_input, 'Protein was not found in list of proteins')
            if not protein.fasta_content:
                raise NoLabsException(ErrorCodes.protein_fasta_is_empty, 'Cannot run folding on a protein without fasta content')

        self.backend = backend
        self.proteins = proteins

    def result_valid(self) -> bool:
        return not not self.foldings

    def set_result(self, result: List[Tuple[Protein, str | bytes]]):
        if not self.proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if not result:
            raise NoLabsException(ErrorCodes.invalid_job_result)

        self.foldings = []

        for protein, pdb in result:
            if not [p for p in self.proteins if p.iid == protein.iid]:
                raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

            self.foldings.append(
                FoldingJobResult(
                    protein_id=protein.iid.value,
                    pdb_content=pdb if isinstance(pdb, bytes) else pdb.encode('utf-8')
                )
            )
