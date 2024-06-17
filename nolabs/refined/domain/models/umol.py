__all__ = [
    'UmolBindingJob'
]

from typing import List
from uuid import UUID

from mongoengine import ReferenceField, ListField, PULL, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, CASCADE, IntField, BinaryField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Job, Protein, Ligand


class UmolBindingJob(Job):
    # region Inputs

    protein: Protein = ReferenceField(Protein, required=True, reverse_delete_rule=CASCADE)
    ligand: Ligand = ReferenceField(Ligand, required=True, reverse_delete_rule=CASCADE)
    pocket_ids: List[int] = ListField(IntField(), required=True)

    # endregion

    predicted_pdb: bytes = BinaryField(required=False)
    predicted_sdf: bytes = BinaryField(required=False)
    plddt_array: List[int] = ListField(IntField())

    def set_input(self,
                  protein: Protein,
                  ligand: Ligand,
                  pocket_ids: List[int]
                  ):
        if not protein:
            raise NoLabsException(ErrorCodes.protein_is_undefined)

        if not protein.fasta_content:
            raise NoLabsException(ErrorCodes.protein_fasta_is_empty, 'Cannot run umol binding job on a protein that does not have a fasta content specificed')

        if not protein.msa:
            raise NoLabsException(ErrorCodes.protein_msa_is_empty, 'Cannot run umol binding job on a protein that does not have a MSA content specificed')

        if not ligand:
            raise NoLabsException(ErrorCodes.ligand_is_undefined, 'Empty ligands input')

        if not ligand.smiles_content:
            raise NoLabsException(ErrorCodes.ligand_smiles_is_empty, 'Ligand smiles is empty')

        self.protein = protein
        self.ligand = ligand
        self.pocket_ids = pocket_ids

    def result_valid(self) -> bool:
        return not not (self.predicted_pdb or self.predicted_sdf or self.plddt_array)

    def set_result(self,
                   protein: Protein,
                   ligand: Ligand,
                   predicted_pdb: bytes | str,
                   predicted_sdf: bytes | str,
                   plddt_array: List[int]):
        if not protein:
            raise NoLabsException(ErrorCodes.protein_is_undefined)

        if protein != self.protein:
            raise NoLabsException(ErrorCodes.protein_not_found_in_job_inputs)

        if not ligand:
            raise NoLabsException(ErrorCodes.ligand_is_undefined)

        if ligand != self.ligand:
            raise NoLabsException(ErrorCodes.ligand_not_found_in_job_inputs)

        self.predicted_pdb = predicted_pdb.decode('utf-8') if isinstance(predicted_pdb, bytes) else predicted_pdb
        self.predicted_sdf = predicted_sdf.decode('utf-8') if isinstance(predicted_sdf, bytes) else predicted_sdf
        self.plddt_array = plddt_array

