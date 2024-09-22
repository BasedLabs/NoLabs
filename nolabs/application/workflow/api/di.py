from __future__ import annotations

__all__ = ["WorkflowDependencies"]

from nolabs.workflow.api.use_cases import (AllWorkflowSchemasFeature,
                                           CreateWorkflowSchemaFeature,
                                           DeleteWorkflowSchemaFeature,
                                           GetComponentStateFeature,
                                           GetJobStateFeature,
                                           GetWorkflowSchemaFeature,
                                           ResetWorkflowFeature,
                                           StartWorkflowComponentFeature,
                                           StartWorkflowFeature,
                                           UpdateWorkflowSchemaFeature)


class WorkflowDependencies:
    @staticmethod
    def create_workflow_schema() -> CreateWorkflowSchemaFeature:
        return CreateWorkflowSchemaFeature()

    @staticmethod
    def delete_workflow_schema() -> DeleteWorkflowSchemaFeature:
        return DeleteWorkflowSchemaFeature()

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
    def reset_workflow():
        return ResetWorkflowFeature()

    @staticmethod
    def get_job_state():
        return GetJobStateFeature()

    @staticmethod
    def all_workflow_schemas():
        return AllWorkflowSchemasFeature()

    @staticmethod
    def get_component_state() -> GetComponentStateFeature:
        return GetComponentStateFeature()
