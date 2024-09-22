import asyncio
import uuid
from pathlib import Path
from typing import List, Optional, Type

from domain.exceptions import ErrorCodes, NoLabsException
from domain.models.common import Experiment
from microservices.reinvent.service.api_models import \
    RunReinforcementLearningRequest
from prefect import State
from prefect.client.schemas.objects import R
from prefect.states import Cancelled, Completed, Failed
from pydantic import BaseModel

from nolabs.application.small_molecules_design.services import (
    ReinventParametersSaver, ReinventSmilesRetriever)
from nolabs.domain.models.common import (DesignedLigandScore,
                                         DrugLikenessScore, JobId, JobName,
                                         Ligand, LigandName, Protein)
from nolabs.domain.models.small_molecules_design import SmallMoleculesDesignJob
from nolabs.infrastructure.cel import cel as celery
from nolabs.infrastructure.settings import settings
from application.workflow import ComponentFlow
from application.workflow.component import Component


class SmallMoleculesDesignLearningInput(BaseModel):
    proteins_with_pdb: List[uuid.UUID]


class SmallMoleculesDesignLearningOutputItem(BaseModel):
    protein: uuid.UUID
    ligands: List[uuid.UUID]


class SmallMoleculesDesignLearningOutput(BaseModel):
    protein_ligands_pairs: List[SmallMoleculesDesignLearningOutputItem]


class SmallMoleculesDesignFlow(ComponentFlow):
    async def get_jobs(self, inp: SmallMoleculesDesignLearningInput) -> List[uuid.UUID]:
        experiment = Experiment.objects.with_id(self.experiment_id)

        job_ids = []

        for protein_id in inp.proteins_with_pdb:
            protein = Protein.objects.with_id(protein_id)

            job_id = JobId(uuid.uuid4())

            job = SmallMoleculesDesignJob(
                id=job_id,
                name=JobName(f"Small molecules design for {protein.name}"),
                experiment=experiment,
            )

            if not protein.pdb_content:
                raise NoLabsException(
                    ErrorCodes.protein_pdb_is_empty,
                    message="Protein pdb content is undefined",
                )

            job.set_protein(protein=protein)

            job.set_inputs(
                protein=protein,
                center_x=0,
                center_y=0,
                center_z=0,
                size_x=5.0,
                size_y=5.0,
                size_z=5.0,
                batch_size=50,
                minscore=0.4,
                epochs=128,
                throw=False,
            )

            job_dir: Path = settings.reinvent_directory / str(job.id)

            tmp_file = job_dir / (str(job.id) + "tmp.pdb")
            tmp_file.write_bytes(job.protein.pdb_content)

            parameters_saver = ReinventParametersSaver()
            await parameters_saver.save_params(
                job_dir=job_dir, job=job, pdb=tmp_file.read_bytes()
            )

            job.change_sampling_size(5)

            await job.save(cascade=True)

            job_ids.append(job.id)

        return job_ids

    async def job_task(self, job_id: uuid.UUID) -> State[R]:
        job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(job_id)

        if not job:
            return Cancelled(message="Job was not found")

        input_errors = job.input_errors(throw=False)

        if input_errors:
            message = ", ".join(i.message for i in input_errors)
            return Cancelled(message=message)

        if not job.celery_task_id:
            (res, task_id) = celery.reinvent_run_learning(
                request=RunReinforcementLearningRequest(config_id=str(job_id)),
                wait=False,
            )
            job.set_task_id(task_id=task_id)
            await job.save()

        if job.celery_task_id:
            celery_task_id = job.celery_task_id
            result = celery.task_result(celery_task_id)

            while not result.failed() and not result.ready():
                await asyncio.sleep(10.0)
                result = celery.task_result(celery_task_id)

            job.celery_task_id = None
            await job.save()

            if result.failed():
                return Failed()

        def normalize_floats(f: float) -> str:
            return str(f).replace(".", "dot")

        smiles_retriever = ReinventSmilesRetriever()

        for smi in smiles_retriever.get_ligands(job.id):
            name = (
                str(job.protein.name)
                + "-binder-drug-likeness-"
                + normalize_floats(smi.drug_likeness)
                + "-score-"
                + normalize_floats(smi.score)
                + "-stage-"
                + smi.stage
            )
            ligand = Ligand.create(
                experiment=job.experiment,
                name=LigandName(name),
                smiles_content=smi.smiles,
            )
            ligand.set_designed_ligand_score(DesignedLigandScore(smi.score))
            ligand.set_drug_likeness_score(DrugLikenessScore(smi.drug_likeness))
            ligand.save()
            complex = ligand.add_binding(protein=job.protein, confidence=smi.score)
            complex.save()

            job.ligands.append(ligand)

            await job.save()

        return Completed()

    async def gather_jobs(
        self, inp: SmallMoleculesDesignLearningInput, job_ids: List[uuid.UUID]
    ) -> Optional[SmallMoleculesDesignLearningOutput]:
        protein_ligands_pairs = []

        for job_id in job_ids:
            job: SmallMoleculesDesignJob = SmallMoleculesDesignJob.objects.with_id(
                job_id
            )

            if not job:
                self.logger.warning("Could not find job", extra={"job_id": job.id})
                continue

            protein_ligands_pairs.append(
                SmallMoleculesDesignLearningOutputItem(
                    protein=job.protein.id,
                    ligands=[l.id for l in job.ligands],
                )
            )

        return SmallMoleculesDesignLearningOutput(
            protein_ligands_pairs=protein_ligands_pairs
        )


class SmallMoleculesDesignLearningComponent(
    Component[SmallMoleculesDesignLearningInput, SmallMoleculesDesignLearningOutput]
):
    name = "Small molecules design"

    @property
    def component_flow_type(self) -> Type["ComponentFlow"]:
        return SmallMoleculesDesignFlow

    @property
    def input_parameter_type(self) -> Type[SmallMoleculesDesignLearningInput]:
        return SmallMoleculesDesignLearningInput

    @property
    def output_parameter_type(self) -> Type[SmallMoleculesDesignLearningOutput]:
        return SmallMoleculesDesignLearningOutput
