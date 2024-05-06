__all__ = [
    'FoldingJobResult',
    'FoldingJob'
]

from enum import Enum
from typing import List, Tuple
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, BinaryField, StringField, EnumField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, LocalisationProbability


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
    proteins: List[Protein] = ListField(ReferenceField(Protein, required=False, reverse_delete_rule=PULL))
    foldings: List[FoldingJobResult] = EmbeddedDocumentListField(FoldingJobResult)
    backend: FoldingBackendEnum = EnumField(FoldingBackendEnum, required=False)

    def set_inputs(self, proteins: List[Protein], backend: FoldingBackendEnum):
        self.foldings = []

        if not proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        for protein in proteins:
            if not Protein.objects.with_id(protein.iid.value):
                raise NoLabsException(ErrorCodes.protein_not_found)

        self.backend = backend
        self.proteins = proteins

    def set_result(self, result: List[Tuple[Protein, str | bytes]]):
        if not self.proteins:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if not result:
            raise NoLabsException(ErrorCodes.invalid_job_result)

        self.foldings = []

        for protein, pdb in result:
            if not [p for p in self.proteins if p.iid == protein.iid]:
                continue

            self.foldings.append(
                FoldingJobResult(
                    protein_id=protein.iid.value,
                    pdb_content=pdb if isinstance(pdb, bytes) else pdb.encode('utf-8')
                )
            )
