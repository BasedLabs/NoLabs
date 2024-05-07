from typing import Annotated

from fastapi import Depends
from reinvent_microservice import ReinventApi

from nolabs.refined.infrastructure.di import InfrastructureDependencies


class SmallMoleculesDesignDependencies:
    @staticmethod
    def delete_job(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> DeleteJobFeature:
        return DeleteJobFeature(api=api)

    @staticmethod
    def get_job_status(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> GetJobStatus:
        return GetJobStatus(api=api)

    @staticmethod
    def get_job(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> GetJobFeature:
        return GetJobFeature(api=api)

    @staticmethod
    def get_job_logs(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> GetJobLogsFeature:
        return GetJobLogsFeature(api=api)

    @staticmethod
    def get_job_smiles(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> GetJobSmilesFeature:
        return GetJobSmilesFeature(api=api)

    @staticmethod
    def setup_job(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> SetupJobFeature:
        return SetupJobFeature(api=api)

    @staticmethod
    def run_learning(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> RunLearningStageJobFeature:
        return RunLearningStageJobFeature(api=api)

    @staticmethod
    def run_sampling(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> RunSamplingStageJobFeature:
        return RunSamplingStageJobFeature(api=api)

    @staticmethod
    def stop_job(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> StopJobFeature:
        return StopJobFeature(api=api)