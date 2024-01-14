from typing import Dict, List

import pydantic


@pydantic.dataclasses.dataclass
class OboNode:
    name: str
    namespace: str
    edges: Dict[str, List[str]]
