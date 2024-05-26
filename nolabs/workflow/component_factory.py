import uuid
from typing import Dict, Optional, Type, List
from uuid import UUID

from nolabs.workflow.component import PythonComponent
from nolabs.workflow.exceptions import WorkflowException


class PythonComponentFactory:
    components: Dict[str, Type[PythonComponent]] = {}

    def __init__(self, components: Dict[str, Type[PythonComponent]]):
        self.components = components

    def create_component(self, name, id: Optional[UUID] = None) -> PythonComponent:
        if name not in self.components:
            raise WorkflowException('Component was not found in components')

        if id:
            return self.components[name](id)

        return self.components[name](uuid.uuid4())

    def all_components(self) -> List[PythonComponent]:
        return list(comp(uuid.uuid4()) for comp in self.components.values())
