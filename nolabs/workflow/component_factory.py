import uuid
from typing import Dict, Optional, Type, List
from uuid import UUID

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.workflow.component import PythonComponent


class PythonComponentFactory:
    components: Dict[str, Type[PythonComponent]] = {}

    def __init__(self, components: Dict[str, Type[PythonComponent]]):
        self.components = components

    def component_exists(self, name: str) -> bool:
        if name not in self.components:
            return False

        return True

    def create_component_instance(self, name: str, id: Optional[UUID] = None) -> PythonComponent:
        if not self.component_exists(name=name):
            raise NoLabsException(ErrorCodes.component_not_found)

        if id:
            return self.components[name](id)

        return self.components[name](uuid.uuid4())

    def all_components(self) -> List[PythonComponent]:
        return list(comp(uuid.uuid4()) for comp in self.components.values())
