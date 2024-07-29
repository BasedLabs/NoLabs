__all__ = [
    'BlastJob'
]

import datetime
from typing import List
from uuid import UUID

from mongoengine import ReferenceField, EmbeddedDocument, FloatField, EmbeddedDocumentListField, \
    UUIDField, CASCADE, IntField, StringField

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.common import Job, Protein, JobInputError


class Hsp(EmbeddedDocument):
    num = IntField(required=True)
    bit_score = FloatField(required=True)
    score = IntField(required=True)
    evalue = FloatField(required=True)
    query_from = IntField(required=True)
    query_to = IntField(required=True)
    hit_from = IntField(required=True)
    hit_to = IntField(required=True)
    query_frame = IntField(required=True)
    hit_frame = IntField(required=True)
    identity = IntField(required=True)
    positive = IntField(required=True)
    gaps = IntField(required=True)
    align_len = IntField(required=True)
    qseq = StringField(required=True)
    hseq = StringField(required=True)
    midline = StringField(required=True)


class Hit(EmbeddedDocument):
    num = IntField(required=True)
    id = StringField(required=True)
    definition = StringField(required=True)
    accession = StringField(required=True)
    length = IntField(required=True)
    hsps = EmbeddedDocumentListField(Hsp, required=True)


class BlastJobResult(EmbeddedDocument):
    protein_id = UUIDField(required=True)
    program = StringField(required=True)
    database = StringField(required=True)
    query_id = StringField(required=True)
    query_def = StringField(required=True)
    query_len = IntField(required=True)
    hits = EmbeddedDocumentListField(Hit, required=True)

    @staticmethod
    def create(protein_id: UUID, program: str, database: str, query_id: str, query_def: str, query_len: int,
               hits: List[Hit]):
        return BlastJobResult(
            protein_id=protein_id,
            program=program,
            database=database,
            query_id=query_id,
            query_def=query_def,
            query_len=query_len,
            hits=hits
        )


class BlastJob(Job):
    # region Inputs

    protein: Protein = ReferenceField(Protein, reverse_delete_rule=CASCADE, required=True)
    job_type = StringField(required=True)
    descriptions = IntField(default=10)
    alignments = IntField(default=10)
    hitlist_size = IntField(default=10)
    expect = FloatField(default=10.0)

    # endregion

    result: List[BlastJobResult] = EmbeddedDocumentListField(BlastJobResult)

    def set_input(self,
                  protein: Protein,
                  job_type: str,
                  descriptions: int = 10,
                  alignments: int = 10,
                  hitlist_size: int = 10,
                  expect: float = 10.0):
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
        return bool(self.protein and self.job_type)

    def set_result(self, result: List[BlastJobResult]):
        if not self.protein:
            raise NoLabsException(ErrorCodes.invalid_job_input, 'Cannot set a result on a job without inputs')

        if not result:
            raise NoLabsException(ErrorCodes.invalid_job_result, 'Result is empty')

        self.result = result

    def _input_errors(self) -> List[JobInputError]:
        errors = []

        if not self.protein:
            errors.append(
                JobInputError(
                    message='Protein is undefined',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        if self.protein and not self.protein.fasta_content:
            errors.append(
                JobInputError(
                    message='Protein fasta content is undefined',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        if not self.job_type:
            errors.append(
                JobInputError(
                    message='BLAST type is undefined',
                    error_code=ErrorCodes.invalid_job_input
                )
            )

        return errors


# Function to create hits and hsps from BLAST output
def create_hits_and_hsps(blast_output: dict) -> List[Hit]:
    hits = []
    iteration_hits = blast_output['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']
    for hit in iteration_hits:
        hsps = []
        for hsp in hit['Hit_hsps']['Hsp']:
            hsps.append(Hsp(
                num=hsp['Hsp_num'],
                bit_score=hsp['Hsp_bit-score'],
                score=hsp['Hsp_score'],
                evalue=hsp['Hsp_evalue'],
                query_from=hsp['Hsp_query-from'],
                query_to=hsp['Hsp_query-to'],
                hit_from=hsp['Hsp_hit-from'],
                hit_to=hsp['Hsp_hit-to'],
                query_frame=hsp['Hsp_query-frame'],
                hit_frame=hsp['Hsp_hit-frame'],
                identity=hsp['Hsp_identity'],
                positive=hsp['Hsp_positive'],
                gaps=hsp['Hsp_gaps'],
                align_len=hsp['Hsp_align-len'],
                qseq=hsp['Hsp_qseq'],
                hseq=hsp['Hsp_hseq'],
                midline=hsp['Hsp_midline']
            ))
        hits.append(Hit(
            num=hit['Hit_num'],
            id=hit['Hit_id'],
            definition=hit['Hit_def'],
            accession=hit['Hit_accession'],
            length=hit['Hit_len'],
            hsps=hsps
        ))
    return hits
