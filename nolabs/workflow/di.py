__all__ = [
    'WorkflowDependencies'
]

from typing import List

from nolabs.workflow.component import Component


class WorkflowDependencies:
    @staticmethod
    def set_workflow_dependencies() -> PlusOneComponent:
        return PlusOneComponent()

