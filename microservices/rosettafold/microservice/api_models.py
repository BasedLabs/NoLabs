from __future__ import annotations

from typing import List

import pydantic

from microservice.mixins import BaseModelMixin

@pydantic.dataclasses.dataclass
class RunRosettaFoldResponse(BaseModelMixin):
    pdb_content: str | None
    errors: List[str]


@pydantic.dataclasses.dataclass
class IsJobRunningResponse:
    is_running: bool
