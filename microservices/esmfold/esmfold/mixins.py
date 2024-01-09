import dataclasses
from typing import Dict, Any, List


class BaseModelMixin:
    def as_log_dict(self) -> Dict[str, Any]:
        d = self.__dict__
        return d


@dataclasses.dataclass
class ErrorResponseMixing:
    errors: List[str] = dataclasses.field(default_factory=list)