__all__ = ["GetJobFeature", "RunJobFeature", "SetupJobFeature"]

from typing import List
from uuid import UUID, uuid4

from nolabs.application.proteinmpnn.worker_models import (
    RunProteinMPNNPredictionInput,
    RunProteinMPNNPredictionOutput,
)
from nolabs.infrastructure.celery_app_factory import get_celery_app, wait_for_task
from nolabs.infrastructure.mongo_connector import get_connection

from nolabs.application.proteinmpnn.api_models import (
    JobResponse,
    JobResult,
    SetupJobRequest,
)
from nolabs.domain.exceptions import ErrorCodes, NoLabsException
from nolabs.domain.models.common import Experiment, JobId, JobName, Protein
from nolabs.domain.models.proteinmpnn import ProteinMPNNJob, ProteinMPNNResult


def map_job_to_response(job: ProteinMPNNJob) -> JobResponse:
    return JobResponse(
        job_id=job.iid.value,
        job_name=job.name.value,
        num_seq_per_target=job.num_seq_per_target,
        sampling_temp=job.sampling_temp,
        seed=job.seed,
        batch_size=job.batch_size,
        is_homomer=job.is_homomer,
        chains_to_design=job.chains_to_design,
        fixed_positions=job.fixed_positions,
        protein_id=job.protein.iid.value,
        result=[
            JobResult(
                sequence_id=res.id,
                sequence=res.sequence,
                fasta_content=res.fasta_content,
                score=res.score,
                global_score=res.global_score,
                T=res.T,
                sample=res.sample,
                seq_recovery=res.seq_recovery,
            )
            for res in job.results
        ],
    )


class GetJobFeature:
    """
    Use case - get job information.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        job_id = JobId(job_id)
        job: ProteinMPNNJob = ProteinMPNNJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        return map_job_to_response(job)


class SetupJobFeature:
    """
    Use case - create new or update existing job.
    """

    async def handle(self, request: SetupJobRequest) -> JobResponse:
        assert request

        job_id = JobId(request.job_id if request.job_id else uuid4())
        job_name = JobName(
            request.job_name if request.job_name else "New ProteinMPNN design job"
        )

        experiment = Experiment.objects.with_id(request.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        job: ProteinMPNNJob = ProteinMPNNJob.objects.with_id(job_id.value)

        if not job:
            job = ProteinMPNNJob(id=job_id, name=job_name, experiment=experiment)

        protein = Protein.objects.with_id(request.protein_id)

        if not protein:
            raise NoLabsException(ErrorCodes.protein_not_found)

        job.set_input(
            protein=protein,
            num_seq_per_target=request.num_seq_per_target,
            sampling_temp=request.sampling_temp,
            seed=request.seed,
            batch_size=request.batch_size,
            is_homomer=request.is_homomer,
            chains_to_design=request.chains_to_design,
            fixed_positions=request.fixed_positions,
        )
        await job.save(cascade=True)

        return map_job_to_response(job)


class RunJobFeature:
    """
    Use case - start job.
    """

    async def handle(self, job_id: UUID) -> JobResponse:
        assert job_id

        job_id = JobId(job_id)
        job: ProteinMPNNJob = ProteinMPNNJob.objects.with_id(job_id.value)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        task_id = uuid4()
        job.set_task_id(task_id=str(task_id))
        await job.save()

        # Prepare the input for the worker
        worker_input = RunProteinMPNNPredictionInput(
            pdb_contents=job.protein.get_pdb(),
            is_homomer=job.is_homomer,
            chains_to_design=job.chains_to_design,
            fixed_positions=job.fixed_positions,
            num_seq_per_target=job.num_seq_per_target,
            sampling_temp=job.sampling_temp,
            seed=job.seed,
            batch_size=job.batch_size,
        )

        # Send task to the worker
        celery_task = get_celery_app().send_task(
            name='run_protmpnn',
            queue='proteinmpnn',
            kwargs={'param': worker_input.model_dump()}
        )

        # Wait for task completion
        job_result_data = await wait_for_task(task_id=celery_task.id)
        job_result = RunProteinMPNNPredictionOutput(**job_result_data)

        # Process the result and save to database
        db = get_connection()
        session = db.client.start_session()
        with session.start_transaction():
            results = []
            for seq_data in job_result.sequences:
                # Create a new sequence object
                sequence_result = ProteinMPNNResult(
                    id=uuid4(),
                    sequence=seq_data.sequence,
                    fasta_content=seq_data.fasta_content,
                    score=seq_data.score,
                    global_score=seq_data.global_score,
                    T=seq_data.T,
                    sample=seq_data.sample,
                    seq_recovery=seq_data.seq_recovery,
                )
                sequence_result.save(session=session)
                results.append(sequence_result)

            job.set_result(results=results)
            await job.save(session=session, cascade=True)

        return map_job_to_response(job)
