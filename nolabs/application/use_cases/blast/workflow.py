import uuid
from typing import List, Optional, Type

from prefect import State
from prefect.client.schemas.objects import R
from prefect.states import Completed
from pydantic import BaseModel

from nolabs.application.use_cases.blast.services import BlastJobRunner
from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.blast import BlastJob
from nolabs.domain.models.common import Experiment, JobId, JobName, Protein
from application.workflow import ComponentFlow
from application.workflow.component import Component


class BlastComponentInput(BaseModel):
    proteins_with_fasta: List[uuid.UUID]
    descriptions: int = 10
    alignments: int = 10
    hitlist_size: int = 10
    expect: int = 10


class BlastComponentOutput(BaseModel):
    hits: List[str]


class BlastComponent(Component[BlastComponentInput, BlastComponentOutput]):
    name = "Blast"
    description = "Finds similar sequences for proteins and nucleotides."

    @property
    def component_flow_type(self) -> Type["ComponentFlow"]:
        return BlastFlow

    @property
    def input_parameter_type(self) -> Type[BlastComponentInput]:
        return BlastComponentInput

    @property
    def output_parameter_type(self) -> Type[BlastComponentOutput]:
        return BlastComponentOutput


class BlastFlow(ComponentFlow):
    async def get_jobs(self, inp: BlastComponentInput) -> List[uuid.UUID]:
        experiment = Experiment.objects.with_id(self.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job_ids = []

        for protein_id in inp.proteins_with_fasta:
            job = BlastJob(
                id=JobId(uuid.uuid4()),
                name=JobName("New BLAST job"),
                experiment=experiment,
            )

            protein = Protein.objects.with_id(protein_id)

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            job.set_input(
                protein=protein,
                job_type="blastp",
                descriptions=inp.descriptions,
                alignments=inp.alignments,
                hitlist_size=inp.hitlist_size,
                expect=inp.expect,
            )

            await job.save(cascade=True)

            job_ids.append(job.id)

        return job_ids

    async def job_task(self, job_id: uuid.UUID) -> State[R]:
        job: BlastJob = BlastJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        job_runner = BlastJobRunner()
        await job_runner.run(job=job)

        return Completed()

    async def gather_jobs(
        self, inp: BlastComponentInput, job_ids: List[uuid.UUID]
    ) -> Optional[BlastComponentOutput]:
        result = []

        for job_id in job_ids:
            job: BlastJob = BlastJob.objects.with_id(job_id)
            for res in job.result:
                for hit in res.hits:
                    result.append(hit.id)

        return BlastComponentOutput(hits=result)
