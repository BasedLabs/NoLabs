from __future__ import annotations

__all__ = [
    'WorkflowDependencies'
]

from typing import Annotated, Type, Dict

from fastapi import Depends

from nolabs.refined.application.use_cases.folding.workflow import FoldingComponent
from nolabs.workflow.application.use_cases import CreateWorkflowSchemaFeature, GetWorkflowSchemaFeature, \
    SetWorkflowSchemaFeature, StartWorkflowFeature, DeleteWorkflowSchemaFeature, AllWorkflowSchemasFeature
from nolabs.workflow.component import PythonComponent


class WorkflowDependencies:
    @staticmethod
    def available_components() -> Dict[str, Type[PythonComponent]]:
        return {
            FoldingComponent.name: FoldingComponent
        }

    @staticmethod
    def create_workflow_schema(components: Annotated[
        Dict[str, Type[PythonComponent]], Depends(WorkflowDependencies.available_components)]) -> CreateWorkflowSchemaFeature:
        return CreateWorkflowSchemaFeature(
            available_components=components
        )

    @staticmethod
    def delete_workflow_schema() -> DeleteWorkflowSchemaFeature:
        return DeleteWorkflowSchemaFeature()

    @staticmethod
    def get_workflow_schema() -> GetWorkflowSchemaFeature:
        return GetWorkflowSchemaFeature()

    @staticmethod
    def set_workflow_schema(components: Annotated[
        Dict[str, Type[PythonComponent]], Depends(WorkflowDependencies.available_components)]) -> SetWorkflowSchemaFeature:
        return SetWorkflowSchemaFeature(
            available_components=components
        )

    @staticmethod
    def start_workflow(components: Annotated[
        Dict[str, Type[PythonComponent]], Depends(WorkflowDependencies.available_components)]) -> StartWorkflowFeature:
        return StartWorkflowFeature(
            available_components=components
        )

    @staticmethod
    def all_workflow_schemas():
        return AllWorkflowSchemasFeature()
