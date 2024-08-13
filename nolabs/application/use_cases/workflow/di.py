from __future__ import annotations

__all__ = [
    'WorkflowDependencies'
]

import logging
from typing import Annotated, Type, Dict

from fastapi import Depends

from nolabs.application.use_cases.test1.workflow import Test1Component
from nolabs.application.use_cases.test2.workflow import Test2Component
from nolabs.application.use_cases.workflow.use_cases import CreateWorkflowDefinitionFeature, \
    GetWorkflowDefinitionFeature, \
    UpdateWorkflowSchemaFeature, StartWorkflowFeature, DeleteWorkflowDefinitionFeature, AllWorkflowDefinitionsFeature, \
    GetComponentStateFeature, ResetWorkflowFeature, StartWorkflowComponentFeature
from nolabs.application.workflow.component import Component
from nolabs.infrastructure.di import InfrastructureDependencies


class WorkflowDependencies:
    @staticmethod
    def available_components() -> Dict[str, Type[Component]]:
        return {
            Test1Component.name: Test1Component,
            Test2Component.name: Test2Component
        }

    @staticmethod
    def create_workflow_schema(components: Annotated[
        Dict[str, Type[Component]], Depends(WorkflowDependencies.available_components)]) -> CreateWorkflowDefinitionFeature:
        return CreateWorkflowDefinitionFeature(
            available_components=components
        )

    @staticmethod
    def delete_workflow_schema() -> DeleteWorkflowDefinitionFeature:
        return DeleteWorkflowDefinitionFeature()

    @staticmethod
    def get_workflow_schema() -> GetWorkflowDefinitionFeature:
        return GetWorkflowDefinitionFeature()

    @staticmethod
    def update_workflow_schema(components: Annotated[
        Dict[str, Type[Component]], Depends(WorkflowDependencies.available_components)]) -> UpdateWorkflowSchemaFeature:
        return UpdateWorkflowSchemaFeature(
            available_components=components,
        )

    @staticmethod
    def start_workflow(components: Annotated[
        Dict[str, Type[Component]], Depends(WorkflowDependencies.available_components)],
                       logger: Annotated[logging.Logger, Depends(InfrastructureDependencies.logger)]) -> StartWorkflowFeature:
        return StartWorkflowFeature(
            available_components=components,
            logger=logger
        )

    @staticmethod
    def start_workflow_component(components: Annotated[
        Dict[str, Type[Component]], Depends(WorkflowDependencies.available_components)],
                                 logger: Annotated[logging.Logger, Depends(InfrastructureDependencies.logger)]) -> StartWorkflowComponentFeature:
        return StartWorkflowComponentFeature(
            available_components=components,
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
