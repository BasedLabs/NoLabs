from typing import Dict, List

from pydantic import BaseModel


class OboNode(BaseModel):
    name: str
    namespace: str
    edges: Dict[str, List[str]]
