from typing import Annotated

from fastapi import Depends
from reinvent_microservice import ReinventApi

from nolabs.refined.application.small_molecules_design.use_cases import *
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
    def save_properties(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> SavePropertiesFeature:
        return SavePropertiesFeature(api=api)

    @staticmethod
    def start_learning(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> StartLearningJobFeature:
        return StartLearningJobFeature(api=api)

    @staticmethod
    def start_sampling(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> StartSamplingJobFeature:
        return StartSamplingJobFeature(api=api)

    @staticmethod
    def stop_job(api: Annotated[ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)]) -> StopJobFeature:
        return StopJobFeature(api=api)