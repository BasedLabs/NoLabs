__all__ = ["BlastJob"]

import datetime
from typing import List
from uuid import UUID

from mongoengine import (
    CASCADE,
    EmbeddedDocument,
    EmbeddedDocumentListField,
    FloatField,
    IntField,
    ReferenceField,
    StringField,
    UUIDField,
)

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Job, JobInputError, Protein


class Hsp(EmbeddedDocument):
    num: int = IntField(required=True)
    bit_score: float = FloatField(required=True)
    score: int = IntField(required=True)
    evalue: float = FloatField(required=True)
    query_from: int = IntField(required=True)
    query_to: int = IntField(required=True)
    hit_from: int = IntField(required=True)
    hit_to: int = IntField(required=True)
    query_frame: int = IntField(required=True)
    hit_frame: int = IntField(required=True)
    identity: int = IntField(required=True)
    positive: int = IntField(required=True)
    gaps: int = IntField(required=True)
    align_len: int = IntField(required=True)
    qseq: str = StringField(required=True)
    hseq: str = StringField(required=True)
    midline: str = StringField(required=True)


class Hit(EmbeddedDocument):
    num: int = IntField(required=True)
    id: str = StringField(required=True)
    definition: str = StringField(required=True)
    accession: str = StringField(required=True)
    length: int = IntField(required=True)
    hsps = EmbeddedDocumentListField(Hsp, required=True)


class BlastJobResult(EmbeddedDocument):
    protein_id: UUID = UUIDField(required=True)
    program: str = StringField(required=True)
    database: str = StringField(required=True)
    query_id: str = StringField(required=True)
    query_def: str = StringField(required=True)
    query_len: int = IntField(required=True)
    hits = EmbeddedDocumentListField(Hit, required=True)

    @staticmethod
    def create(
        protein_id: UUID,
        program: str,
        database: str,
        query_id: str,
        query_def: str,
        query_len: int,
        hits: List[Hit],
    ):
        return BlastJobResult(
            protein_id=protein_id,
            program=program,
            database=database,
            query_id=query_id,
            query_def=query_def,
            query_len=query_len,
            hits=hits,
        )


class BlastJob(Job):
    # region Inputs

    protein: Protein = ReferenceField(
        Protein, reverse_delete_rule=CASCADE, required=True
    )
    job_type = StringField(required=True)
    descriptions = IntField(default=10)
    alignments = IntField(default=10)
    hitlist_size = IntField(default=10)
    expect = FloatField(default=10.0)

    # endregion

    result: List[BlastJobResult] = EmbeddedDocumentListField(BlastJobResult)

    def set_input(
        self,
        protein: Protein,
        job_type: str,
        descriptions: int = 10,
        alignments: int = 10,
        hitlist_size: int = 10,
        expect: float = 10.0,
    ):
        self.result = []

        self.protein = protein
        self.job_type = job_type
        self.descriptions = descriptions
        self.alignments = alignments
        self.hitlist_size = hitlist_size
        self.expect = expect

        self.input_errors(throw=True)

        self.inputs_updated_at = datetime.datetime.utcnow()

    def result_valid(self) -> bool:
        return bool(self.result)

    def input_valid(self) -> bool:
        # TODO: add more blast jobs types when DNA is introduced, move to enums then
        return bool(self.protein and self.job_type and (self.job_type in ["blastp"]))

    def set_result(self, result: List[BlastJobResult]):
        if not self.protein:
            raise NoLabsException(
                ErrorCodes.invalid_job_input,
                "Cannot set a result on a job without inputs",
            )

        if not result:
            raise NoLabsException(ErrorCodes.invalid_job_result, "Result is empty")

        self.result = result

    def _input_errors(self) -> List[JobInputError]:
        errors = []

        if not self.protein:
            errors.append(
                JobInputError(
                    message="Protein is undefined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if self.protein and not self.protein.fasta_content:
            errors.append(
                JobInputError(
                    message="Protein fasta content is undefined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        if not self.job_type:
            errors.append(
                JobInputError(
                    message="BLAST type is undefined",
                    error_code=ErrorCodes.invalid_job_input,
                )
            )

        return errors
