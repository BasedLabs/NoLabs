import uuid
from typing import List

from nolabs.application.adaptyv_bio.api_proxy import AdaptyvBioProteinAffinityCharacterizationApi
from nolabs.application.adaptyv_bio.protein_affinity_characterization.api_models import JobResponse, SetupJobRequest, \
    TargetResponse, EstimatesResponse
from nolabs.domain.exceptions import NoLabsException, ErrorCodes
from nolabs.domain.models.adaptyv_bio.protein_affinity_characterization_job import ProteinAffinityCharacterizationJob
from nolabs.workflow.core.graph import Graph


class StartJobFeature:
    async def handle(self, job_id: uuid.UUID):
        job: ProteinAffinityCharacterizationJob = ProteinAffinityCharacterizationJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        job.input_errors(throw=True)

        experiment_id = job.component.experiment.id
        component_id = job.component.id
        graph = Graph(experiment_id=experiment_id)
        await graph.schedule(component_ids=[component_id])


class GetJobFeature:
    async def handle(self, job_id: uuid.UUID) -> JobResponse:
        job: ProteinAffinityCharacterizationJob = ProteinAffinityCharacterizationJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        return JobResponse(
            job_id=job_id,
            job_name=job.name.value,
            number_of_designs=job.number_of_designs,
            dna_length=job.dna_length,
            replicates=job.replicates,
            report_email=job.report_email,
            target_id=job.target_id,
            cart_total=job.cart_total,
            session_url=job.session_url,
            swissprot_id=job.swissprot_id
        )

class SetupJobFeature:
    async def handle(self, request: SetupJobRequest):
        job: ProteinAffinityCharacterizationJob = ProteinAffinityCharacterizationJob.objects.with_id(request.job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        job.set_input(
            number_of_designs=request.number_of_designs,
            dna_length=request.dna_length,
            replicates=request.replicates,
            report_email=request.report_email,
            target_id=request.target_id,
            cart_total=request.cart_total,
            swissprot_id=request.swissprot_id
        )

class ListTargetsFeature:
    async def handle(self, search_query: str) -> List[TargetResponse]:
        api = AdaptyvBioProteinAffinityCharacterizationApi()
        targets = api.list_targets(query=search_query)
        return [TargetResponse(id=t.id,
                               name=t.name,
                               description=t.description,
                               swissprot_id=t.swissprot_id) for t in targets]

class GetEstimatesFeature:
    async def handle(self, job_id: uuid.UUID) -> EstimatesResponse:
        job: ProteinAffinityCharacterizationJob = ProteinAffinityCharacterizationJob.objects.with_id(job_id)

        if not job:
            raise NoLabsException(ErrorCodes.job_not_found)

        api = AdaptyvBioProteinAffinityCharacterizationApi()
        estimates = api.get_experiment_estimates(
            n_designs=job.number_of_designs,
            avg_length=job.dna_length,
            n_replicates=job.replicates
        )
        return EstimatesResponse(
            total_price=estimates.total_price,
            turnaround_time=estimates.turnaround_time
        )