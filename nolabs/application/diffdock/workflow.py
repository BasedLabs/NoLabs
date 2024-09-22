__all__ = ["DiffDockComponent"]

import uuid
from typing import List, Type

from domain.exceptions import ErrorCodes, NoLabsException
from microservices.diffdock.service.api_models import (
    RunDiffDockPredictionRequest)
from prefect import State
from prefect.client.schemas.objects import R
from prefect.states import Cancelled, Completed, Failed
from pydantic import BaseModel

from nolabs.domain.models.common import JobId, JobName, Ligand, Protein
from nolabs.domain.models.diffdock import DiffDockBindingJob, DiffDockJobResult
from nolabs.infrastructure.cel import cel as celery
from application.workflow import ComponentFlow
from application.workflow.component import Component, TInput, TOutput


class DiffDockComponentInput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]
    ligands: List[uuid.UUID]


class DiffDockComponentOutput(BaseModel):
    protein_complexes: List[uuid.UUID]


class DiffDockComponent(Component[DiffDockComponentInput, DiffDockComponentOutput]):
    name = "DiffDock"
    description = "Prediction of protein-ligand complexes"

    @property
    def input_parameter_type(self) -> Type[TInput]:
        return DiffDockComponentInput

    @property
    def output_parameter_type(self) -> Type[TOutput]:
        return DiffDockComponentOutput

    @property
    def component_flow_type(self) -> Type["ComponentFlow"]:
        return DiffdockComponentFlow


class DiffdockComponentFlow(ComponentFlow):
    async def get_jobs(self, inp: DiffDockComponentInput) -> List[uuid.UUID]:
        job_ids = []

        for protein_id in inp.proteins_with_pdb:
            for ligand_id in inp.ligands:
                protein: Protein = Protein.objects.with_id(protein_id)

                if not protein:
                    raise NoLabsException(
                        ErrorCodes.protein_not_found, data={"protein_id": protein_id}
                    )

                ligand: Ligand = Ligand.objects.with_id(ligand_id)

                if not ligand:
                    raise NoLabsException(
                        ErrorCodes.ligand_not_found, data={"ligand_id": ligand_id}
                    )

                job = DiffDockBindingJob.objects(protein=protein.id, ligand=ligand.id).first()

                if not job:
                    job = DiffDockBindingJob(
                        id=JobId(uuid.uuid4()),
                        name=JobName("Diffdock binding job"),
                        experiment=self.experiment_id,
                    )
                    job.set_input(protein=protein, ligand=ligand)
                    await job.save()

                job_ids.append(job.id)

        return job_ids

    async def job_task(self, job_id: uuid.UUID) -> State[R]:
        job: DiffDockBindingJob = DiffDockBindingJob.objects.with_id(job_id)

        input_errors = job.input_errors(throw=False)

        if input_errors:
            message = ", ".join(i.message for i in input_errors)

            return Cancelled(message=message)

        request = RunDiffDockPredictionRequest(
            pdb_contents=job.protein.get_pdb(), sdf_contents=job.ligand.get_sdf()
        )
        (job_result, task_id) = await celery.diffdock_inference(payload=request)

        ligand = job.ligand

        result = []

        for item in job_result.sdf_results:
            complex = ligand.add_binding(
                protein=job.protein,
                sdf_content=item.sdf_content,
                minimized_affinity=item.minimized_affinity,
                scored_affinity=item.scored_affinity,
                confidence=item.confidence,
                plddt_array=[],
                name=item.sdf_file_name,
            )

            complex.save()

            result.append(
                DiffDockJobResult.create(
                    complex_id=complex.iid.value,
                    sdf_content=item.sdf_content,
                    minimized_affinity=item.minimized_affinity,
                    scored_affinity=item.scored_affinity,
                    confidence=item.confidence,
                )
            )

        job.set_result(result=result)

        await job.save()

        if not job.result_valid():
            return Failed(message="Result of job is invalid")

        return Completed()

    async def gather_jobs(self, inp: DiffDockComponentOutput, job_ids: List[uuid.UUID]):
        ids = []

        for job_id in job_ids:
            job: DiffDockBindingJob = DiffDockBindingJob.objects.with_id(job_id)

            for complex in job.result:
                ids.append(complex.complex_id)

        return DiffDockComponentOutput(protein_complexes=ids)
