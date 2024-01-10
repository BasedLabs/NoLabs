import dataclasses
from typing import Dict, List

import pydantic


@dataclasses.dataclass
@pydantic.dataclass
class OboNode:
    name: str
    namespace: str
    edges: Dict[str, List[str]]
