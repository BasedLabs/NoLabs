from fastapi import APIRouter

router = APIRouter(
    prefix='/api/v1/workflow',
    tags=['Workflow'],

)


@router.get('/{experiment_id}', summary='Get workflow schema')
async def get_schema(
        feature: Annotated[
            RunJobFeature, Depends(GeneOntologyDependencies.run_job)],
        job_id: UUID
) -> JobResponse:
    return await feature.handle(job_id=job_id)