import dataclasses
from typing import Dict, Any, List


class BaseModelMixin:
    def as_log_dict(self) -> Dict[str, Any]:
        max_val_len = 40
        d = self.__dict__.copy()
        for key in d:
            if isinstance(d[key], str) and len(d[key]) > max_val_len:
                d[key] = d[key][:max_val_len]
        return d


@dataclasses.dataclass
class ErrorResponseMixing:
    errors: List[str] = dataclasses.field(default_factory=list)