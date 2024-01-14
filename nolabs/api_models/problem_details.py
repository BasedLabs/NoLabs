from typing import List

import pydantic.dataclasses as pcdataclasses

from nolabs.exceptions import ErrorCodes


@pcdataclasses.dataclass
class ProblemDetailsResponse:
    errors: List[str]
    error_code: ErrorCodes
