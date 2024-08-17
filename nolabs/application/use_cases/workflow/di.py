from __future__ import annotations

__all__ = [
    'WorkflowDependencies'
]

import logging
from typing import Annotated

from fastapi import Depends

from nolabs.application.use_cases.workflow.use_cases import CreateWorkflowDefinitionFeature, \
    GetWorkflowDefinitionFeature, \
    UpdateWorkflowDefinitionFeature, StartWorkflowFeature, DeleteWorkflowDefinitionFeature, \
    AllWorkflowDefinitionsFeature, \
    GetComponentStateFeature, ResetWorkflowFeature, StartWorkflowComponentFeature
from nolabs.infrastructure.di import InfrastructureDependencies


class WorkflowDependencies:
    @staticmethod
    def create_workflow_schema() -> CreateWorkflowDefinitionFeature:
        return CreateWorkflowDefinitionFeature(
        )

    @staticmethod
    def delete_workflow_schema() -> DeleteWorkflowDefinitionFeature:
        return DeleteWorkflowDefinitionFeature()

    @staticmethod
    def get_workflow_schema() -> GetWorkflowDefinitionFeature:
        return GetWorkflowDefinitionFeature()

    @staticmethod
    def update_workflow_schema() -> UpdateWorkflowDefinitionFeature:
        return UpdateWorkflowDefinitionFeature(
        )

    @staticmethod
    def start_workflow() -> StartWorkflowFeature:
        return StartWorkflowFeature()

    @staticmethod
    def start_workflow_component(logger: Annotated[logging.Logger, Depends(InfrastructureDependencies.logger)]) -> StartWorkflowComponentFeature:
        return StartWorkflowComponentFeature(
            logger=logger
        )

    @staticmethod
    def reset_workflow():
        return ResetWorkflowFeature()

    @staticmethod
    def all_workflow_schemas():
        return AllWorkflowDefinitionsFeature()

    @staticmethod
    def get_component_state() -> GetComponentStateFeature:
        return GetComponentStateFeature()
