from __future__ import annotations

__all__ = ["WorkflowDependencies"]

from nolabs.application.workflow.api.use_cases import (
    CreateWorkflowSchemaFeature,
    GetComponentStateFeature,
    GetJobStateFeature,
    GetWorkflowSchemaFeature,
    StartWorkflowComponentFeature,
    StartWorkflowFeature,
    UpdateWorkflowSchemaFeature,
)


class WorkflowDependencies:
    @staticmethod
    def create_workflow_schema() -> CreateWorkflowSchemaFeature:
        return CreateWorkflowSchemaFeature()

    @staticmethod
    def get_workflow_schema() -> GetWorkflowSchemaFeature:
        return GetWorkflowSchemaFeature()

    @staticmethod
    def update_workflow_schema() -> UpdateWorkflowSchemaFeature:
        return UpdateWorkflowSchemaFeature()

    @staticmethod
    def start_workflow() -> StartWorkflowFeature:
        return StartWorkflowFeature()

    @staticmethod
    def start_workflow_component() -> StartWorkflowComponentFeature:
        return StartWorkflowComponentFeature()

    @staticmethod
    def get_job_state():
        return GetJobStateFeature()

    @staticmethod
    def get_component_state() -> GetComponentStateFeature:
        return GetComponentStateFeature()
