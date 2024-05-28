from __future__ import annotations

__all__ = [
    'WorkflowDependencies'
]

from typing import Annotated

from fastapi import Depends

from nolabs.refined.application.use_cases.folding.workflow import FoldingComponent
from nolabs.workflow.application.use_cases import CreateWorkflowSchemaFeature, GetWorkflowSchemaFeature, \
    SetWorkflowSchemaFeature, StartWorkflowFeature, StopWorkflowFeature, GetComponentJobIdsFeature
from nolabs.workflow.component_factory import PythonComponentFactory


class WorkflowDependencies:
    @staticmethod
    def components_factory() -> PythonComponentFactory:
        return PythonComponentFactory(
            components={
                FoldingComponent.name: FoldingComponent
            }
        )

    @staticmethod
    def create_workflow_schema(factory: Annotated[PythonComponentFactory, Depends(WorkflowDependencies.components_factory)]) -> CreateWorkflowSchemaFeature:
        return CreateWorkflowSchemaFeature(
            factory=factory
        )

    @staticmethod
    def get_workflow_schema() -> GetWorkflowSchemaFeature:
        return GetWorkflowSchemaFeature()

    @staticmethod
    def set_workflow_schema(factory: Annotated[PythonComponentFactory, Depends(WorkflowDependencies.components_factory)]) -> SetWorkflowSchemaFeature:
        return SetWorkflowSchemaFeature(
            factory=factory
        )

    @staticmethod
    def start_workflow(factory: Annotated[
        PythonComponentFactory, Depends(WorkflowDependencies.components_factory)]) -> StartWorkflowFeature:
        return StartWorkflowFeature(
            factory=factory
        )

    @staticmethod
    def stop_workflow(factory: Annotated[
        PythonComponentFactory, Depends(WorkflowDependencies.components_factory)]) -> StopWorkflowFeature:
        return StopWorkflowFeature(
            factory=factory
        )

    @staticmethod
    def get_component_job_ids(factory: Annotated[PythonComponentFactory, Depends(WorkflowDependencies.components_factory)]) -> GetComponentJobIdsFeature:
        return GetComponentJobIdsFeature(
            factory=factory
        )
