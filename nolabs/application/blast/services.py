from typing import List

from nolabs.application.blast import blast_api
from nolabs.application.blast.blast_api import BlastType
from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.blast import BlastJob, BlastJobResult, Hit, Hsp


class BlastJobRunner:
    async def run(self, job: BlastJob):
        blast_output = await blast_api.run_blast(
            program=BlastType.BLASTP,
            sequence=job.protein.get_amino_acid_sequence(),
            descriptions=job.descriptions,
            alignments=job.alignments,
            hitlist_size=job.hitlist_size,
            expect=job.expect,
        )

        if not blast_output:
            raise NoLabsException(ErrorCodes.blast_api_error)

        hits = self._create_hits_and_hsps(blast_output)

        blast_job_result = BlastJobResult.create(
            protein_id=job.protein.id,
            program=blast_output["BlastOutput"]["BlastOutput_program"],
            database=blast_output["BlastOutput"]["BlastOutput_db"],
            query_id=blast_output["BlastOutput"]["BlastOutput_query-ID"],
            query_def=blast_output["BlastOutput"]["BlastOutput_query-def"],
            query_len=int(blast_output["BlastOutput"]["BlastOutput_query-len"]),
            hits=hits,
        )

        job.set_result([blast_job_result])

        await job.save()

    def _create_hits_and_hsps(self, blast_output: dict) -> List[Hit]:
        hits = []
        iteration_hits = blast_output["BlastOutput"]["BlastOutput_iterations"][
            "Iteration"
        ]["Iteration_hits"]["Hit"]
        for hit in iteration_hits:
            hsps = []
            hsp = hit["Hit_hsps"]["Hsp"]
            hsps.append(
                Hsp(
                    num=int(hsp["Hsp_num"]),
                    bit_score=float(hsp["Hsp_bit-score"]),
                    score=int(hsp["Hsp_score"]),
                    evalue=float(hsp["Hsp_evalue"]),
                    query_from=int(hsp["Hsp_query-from"]),
                    query_to=int(hsp["Hsp_query-to"]),
                    hit_from=int(hsp["Hsp_hit-from"]),
                    hit_to=int(hsp["Hsp_hit-to"]),
                    query_frame=int(hsp["Hsp_query-frame"]),
                    hit_frame=int(hsp["Hsp_hit-frame"]),
                    identity=int(hsp["Hsp_identity"]),
                    positive=int(hsp["Hsp_positive"]),
                    gaps=int(hsp["Hsp_gaps"]),
                    align_len=int(hsp["Hsp_align-len"]),
                    qseq=hsp["Hsp_qseq"],
                    hseq=hsp["Hsp_hseq"],
                    midline=hsp["Hsp_midline"],
                )
            )
            hits.append(
                Hit(
                    num=int(hit["Hit_num"]),
                    id=hit["Hit_id"],
                    definition=hit["Hit_def"],
                    accession=hit["Hit_accession"],
                    length=int(hit["Hit_len"]),
                    hsps=hsps,
                )
            )
        return hits
