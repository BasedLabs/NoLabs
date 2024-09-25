__all__ = ["GetJobFeature", "RunJobFeature", "SetupJobFeature"]

import uuid
from uuid import UUID

from microservices.esmfold_light.service.api_models import InferenceInput
from mongoengine import Q

from nolabs.application.folding.api_models import (JobResponse, JobResult,
                                                   SetupJobRequest)
from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Experiment, JobId, JobName, Protein
from nolabs.domain.models.folding import FoldingBackendEnum, FoldingJob
from nolabs.infrastructure.cel import cel as celery
from nolabs.infrastructure.log import logger
from nolabs.utils import generate_uuid


def map_job_to_response(job: FoldingJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        backend=FoldingBackendEnum(job.backend),
        protein_id=job.protein.id,
        result=JobResult(
            protein_id=job.folded_protein.id, pdb=job.folded_protein.get_pdb()
        ),
        experiment_id=job.experiment.id,
    )


class GetJobFeature:
    async def handle(self, job_id: UUID) -> JobResponse:
        extra = {
            "job_id": job_id,
        }

        logger.info("Get folding job", extra=extra)

        try:
            job_id = JobId(job_id)
            job: FoldingJob = FoldingJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            logger.info("Get folding job success", extra=extra)

            return map_job_to_response(job)
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e
            raise NoLabsException(ErrorCodes.get_folding_job_failed) from e


class SetupJobFeature:
    async def handle(self, request: SetupJobRequest) -> JobResponse:
        extra = {
            "job_id": request.job_id,
        }

        logger.info("Setup folding job", extra=extra)

        try:
            job_id = JobId(request.job_id if request.job_id else generate_uuid())
            job_name = JobName(
                request.job_name if request.job_name else "New folding job"
            )
            folding_backend = (
                FoldingBackendEnum(request.backend)
                if request.backend
                else FoldingBackendEnum.esmfold
            )

            experiment = Experiment.objects.with_id(request.experiment_id)

            if not experiment:
                raise NoLabsException(ErrorCodes.experiment_not_found)

            job: FoldingJob = FoldingJob.objects(
                Q(id=job_id.value) | Q(name=job_name.value)
            ).first()

            if not job:
                job: FoldingJob = FoldingJob.create(
                    id=job_id, name=job_name, experiment=experiment
                )

            protein = Protein.objects.with_id(request.protein_id)

            if not protein:
                raise NoLabsException(ErrorCodes.protein_not_found)

            job.set_inputs(protein=protein, backend=folding_backend)
            await job.save(cascade=True)

            logger.info("Setup folding job success", extra=extra)

            return map_job_to_response(job)
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e
            raise NoLabsException(ErrorCodes.setup_folding_job_failed) from e


class RunJobFeature:
    async def handle(self, job_id: UUID) -> JobResponse:
        extra = {"job_id": job_id}

        logger.info("Setup folding job", extra=extra)

        try:
            job_id = JobId(job_id)
            job: FoldingJob = FoldingJob.objects.with_id(job_id.value)

            if not job:
                raise NoLabsException(ErrorCodes.job_not_found)

            protein = job.protein

            sequence = protein.get_amino_acid_sequence()

            if not sequence:
                raise NoLabsException(
                    ErrorCodes.protein_fasta_is_empty,
                    f"Protein sequence is empty for protein {protein.name}",
                )

            if job.backend == FoldingBackendEnum.esmfold_light:
                inference_result = await celery.esmfold_light_inference(
                    task_id=job.id,
                    payload=InferenceInput(fasta_sequence=protein.get_fasta()),
                )
                id = uuid.uuid4()
                folded_protein = protein.copy(id=id)
                folded_protein.set_pdb(inference_result.pdb_content)
                job.set_result(protein=folded_protein)
                await job.save(cascade=True)

                return map_job_to_response(job)

            raise NoLabsException(ErrorCodes.folding_method_unknown)
        except Exception as e:
            if isinstance(e, NoLabsException):
                raise e
            raise NoLabsException(ErrorCodes.run_folding_job_failed) from e
