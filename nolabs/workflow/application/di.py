__all__ = [
    'WorkflowDependencies'
]

from typing import List, Annotated

from fastapi import Depends

from nolabs.refined.application.use_cases.test_app.components import PlusOneComponent, PlusTwoComponent
from nolabs.refined.application.use_cases.test_app.di import TestAppDependencies
from nolabs.workflow.function import Component


class WorkflowDependencies:
    @staticmethod
    def components(
            one_plus_one: Annotated[PlusOneComponent, Depends(TestAppDependencies.plus_one)],
            one_plus_two: Annotated[PlusTwoComponent, Depends(TestAppDependencies.plus_two)]
    ) -> List[Component]:
        return [
            one_plus_one,
            one_plus_two
        ]

