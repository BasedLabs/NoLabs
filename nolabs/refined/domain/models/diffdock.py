__all__ = [
    'DiffDockBindingJob'
]

from typing import List
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, CASCADE, IntField, BinaryField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, Ligand


class DiffDockJobResult(EmbeddedDocument):
    ligand_id: UUID = UUIDField(required=True)
    sdf_content: bytes = BinaryField(required=True)
    minimized_affinity: float = FloatField(required=False)
    scored_affinity: float = FloatField(required=False)
    confidence: float = FloatField(required=False)

    def __init__(self,
                 ligand_id: UUID,
                 sdf_content: bytes | str,
                 minimized_affinity: float,
                 scored_affinity: float,
                 confidence: float,
                 *args, **kwargs):
        if isinstance(sdf_content, str):
            sdf_content = sdf_content.encode('utf-8')

        super().__init__(
            ligand_id=ligand_id,
            sdf_content=sdf_content,
            minimized_affinity=minimized_affinity,
            scored_affinity=scored_affinity,
            confidence=confidence,
            *args, **kwargs)


class DiffDockBindingJob(Job):
    protein: Protein = ReferenceField(Protein, required=False, reverse_delete_rule=CASCADE)
    ligands: List[Ligand] = ListField(ReferenceField(Ligand, required=False, reverse_delete_rule=PULL))
    samples_per_complex: int = IntField(required=False, default=40)

    result: List[DiffDockJobResult] = EmbeddedDocumentListField(DiffDockJobResult)

    def set_input(self,
                  protein: Protein,
                  ligands: List[Ligand],
                  samples_per_complex: int
                  ):
        if not protein:
            raise NoLabsException(ErrorCodes.protein_is_undefined)

        if not protein.pdb_content:
            raise NoLabsException(ErrorCodes.protein_pdb_is_empty, 'Cannot run binding job on empty pdb')

        if not ligands:
            raise NoLabsException(ErrorCodes.ligand_is_undefined, 'Empty ligands input')

        if samples_per_complex <= 0:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Samples per complex must be greater than 0')

        self.protein = protein
        self.ligands = ligands
        self.samples_per_complex = samples_per_complex

    def set_result(self,
                   protein: Protein,
                   ligands: List[Ligand],
                   result: List[DiffDockJobResult]):
        if not protein:
            raise NoLabsException(ErrorCodes.invalid_job_input)

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if set(l.iid for l in self.ligands) != set(l.iid for l in ligands):
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Some or all of ligands were not a part of job input')

        if not result:
            raise NoLabsException(ErrorCodes.invalid_job_result, 'Result is empty')

        self.result = result
