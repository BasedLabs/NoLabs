from typing import Annotated

from fastapi import Depends
from reinvent_microservice import ReinventApi

from nolabs.application.use_cases.small_molecules_design.use_cases import (
    DeleteJobFeature, GetJobFeature, GetJobLogsFeature, GetJobSmilesFeature,
    GetJobStatusFeature, RunLearningStageJobFeature,
    RunSamplingStageJobFeature, SetupJobFeature, StopJobFeature)
from nolabs.infrastructure.di import InfrastructureDependencies


class SmallMoleculesDesignDependencies:
    @staticmethod
    def delete_job(
        api: Annotated[
            ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)
        ]
    ) -> DeleteJobFeature:
        return DeleteJobFeature(api=api)

    @staticmethod
    def get_job_status(
        api: Annotated[
            ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)
        ]
    ) -> GetJobStatusFeature:
        return GetJobStatusFeature(api=api)

    @staticmethod
    def get_job(
        api: Annotated[
            ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)
        ]
    ) -> GetJobFeature:
        return GetJobFeature(api=api)

    @staticmethod
    def get_job_logs(
        api: Annotated[
            ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)
        ]
    ) -> GetJobLogsFeature:
        return GetJobLogsFeature(api=api)

    @staticmethod
    def get_job_smiles(
        api: Annotated[
            ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)
        ]
    ) -> GetJobSmilesFeature:
        return GetJobSmilesFeature(api=api)

    @staticmethod
    def setup_job(
        api: Annotated[
            ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)
        ]
    ) -> SetupJobFeature:
        return SetupJobFeature(api=api)

    @staticmethod
    def run_learning(
        api: Annotated[
            ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)
        ]
    ) -> RunLearningStageJobFeature:
        return RunLearningStageJobFeature(api=api)

    @staticmethod
    def run_sampling(
        api: Annotated[
            ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)
        ]
    ) -> RunSamplingStageJobFeature:
        return RunSamplingStageJobFeature(api=api)

    @staticmethod
    def stop_job(
        api: Annotated[
            ReinventApi, Depends(InfrastructureDependencies.reinvent_microservice)
        ]
    ) -> StopJobFeature:
        return StopJobFeature(api=api)
