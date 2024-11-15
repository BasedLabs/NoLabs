__all__ = ["DiffDockComponent"]

import uuid
from typing import List, Type, Optional, Dict, Any

from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from pydantic import BaseModel

from nolabs.application.diffdock.worker_models import RunDiffDockPredictionInput, RunDiffDockPredictionOutput
from nolabs.domain.models.common import JobId, JobName, Ligand, Protein
from nolabs.domain.models.diffdock import DiffDockBindingJob
from nolabs.infrastructure.mongo_connector import get_connection
from nolabs.workflow.core.component import Component
from nolabs.workflow.core.flow import ComponentFlowHandler


class DiffDockComponentInput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]
    ligands: List[uuid.UUID]
    number_of_complexes: int = 3


class DiffDockComponentOutput(BaseModel):
    protein_complexes: List[uuid.UUID]


class DiffDockComponent(Component[DiffDockComponentInput, DiffDockComponentOutput]):
    name = "DiffDock"
    description = "Prediction of protein-ligand complexes"

    @property
    def input_parameter_type(self) -> Type[DiffDockComponentInput]:
        return DiffDockComponentInput

    @property
    def output_parameter_type(self) -> Type[DiffDockComponentOutput]:
        return DiffDockComponentOutput

    @property
    def component_flow_type(self) -> Type["ComponentFlowHandler"]:
        return DiffdockComponentFlowHandler


class DiffdockComponentFlowHandler(ComponentFlowHandler):
    async def on_component_task(self, inp: DiffDockComponentInput) -> List[uuid.UUID]:
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

                job = DiffDockBindingJob.objects(
                    protein=protein.id, ligand=ligand.id
                ).first()

                if not job:
                    job = DiffDockBindingJob.create(
                        id=JobId(uuid.uuid4()),
                        name=JobName("Diffdock binding job"),
                        component=self.component_id
                    )
                    job.set_input(protein=protein, ligand=ligand, samples_per_complex=inp.number_of_complexes)
                    await job.save()

                job_ids.append(job.id)

        return job_ids

    async def on_job_task(self, job_id: uuid.UUID):
        job: DiffDockBindingJob = DiffDockBindingJob.objects.with_id(job_id)

        input_errors = job.input_errors(throw=False)

        if input_errors:
            message = ", ".join(i.message for i in input_errors)

            raise NoLabsException(ErrorCodes.invalid_job_input, message=message)

        request = RunDiffDockPredictionInput(
            pdb_contents=job.protein.get_pdb(),
            sdf_contents=job.ligand.get_sdf(),
            samples_per_complex=job.samples_per_complex
        )
        task_id = uuid.uuid4()
        job.set_task_id(task_id=str(task_id))
        await job.save()

        return await self.schedule(
            job_id=job.id,
            celery_task_name="inference",
            celery_queue="diffdock",
            input={"param": request.model_dump()}
        )

    async def on_job_completion(
        self, job_id: uuid.UUID, long_running_output: Optional[Dict[str, Any]]
    ):
        job: DiffDockBindingJob = DiffDockBindingJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        job_result = RunDiffDockPredictionOutput(**long_running_output)

        ligand = job.ligand
        result = []

        db = get_connection()
        session = db.client.start_session()
        with session.start_transaction():
            for item in job_result.sdf_results:
                ligand_for_complex = Ligand.copy(ligand)

                ligand_for_complex.save(session=session)

                complex = Protein.create_complex(
                    protein=job.protein,
                    ligand=ligand_for_complex,
                    minimized_affinity=item.minimized_affinity,
                    scored_affinity=item.scored_affinity,
                    confidence=item.confidence,
                    plddt_array=[]
                )

                complex.save(session=session)
                result.append(complex)

            job.set_result(complexes=result)

            await job.save(session=session, cascade=True)
