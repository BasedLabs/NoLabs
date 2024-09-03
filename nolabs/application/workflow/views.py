import uuid
from typing import Any, Dict, List

from pydantic import BaseModel


class InputPropertyErrorModel(EmbeddedDocument):
    loc: List[str] = ListField(StringField())
    msg: str = StringField()

    @classmethod
    def create(cls, loc: List[str], msg: str) -> 'InputPropertyErrorModel':
        return InputPropertyErrorModel(
            loc=loc,
            msg=msg
        )


class ComponentRunView(BaseModel):
    input_dict: Dict[str, Any]
    output_dict: Dict[str, Any]
    job_ids: List[uuid.UUID]
    input_property_errors = [],
    last_exceptions = [],
    jobs_errors = []