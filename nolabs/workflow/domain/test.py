import json
from enum import Enum
from typing import Optional, List, Tuple, Dict, Any
from uuid import UUID, uuid4

from fastapi import UploadFile
from pydantic import BaseModel, StrictBytes
from pydantic.dataclasses import dataclass

from nolabs.workflow.schema import Schema


class En2(Enum):
    c = 'c'
class En(Enum):
    a = 'a'
    b = 'b'

class AnotherInput2(BaseModel):
    """Hello there"""
    a: UUID
    vv: int
    f: UploadFile
    b: bytes
    en: En = En.a
    d: Dict[int, int]


class AnotherInput(BaseModel):
    another_value: str
    inp: AnotherInput2
    value: Optional[float] = 10.0


class Input(BaseModel):
    """
    Sum
    :param a: Str input
    :param b: Int b
    """
    a: str
    en2: En2
    b: int
    c: List[Tuple[int, str]]
    parameter: AnotherInput | None = None
    optional: Optional[int] = 10
    d: Dict[str, Any]



#class Output(BaseModel):
#    c: float
#
#input_schema = Schema.get_schema(Input)
#another_input_schema = Schema.get_schema(AnotherInput2)
#
#o = Output(
#    c=10.0
#)
#print(o.dict())

pretty_schema = json.dumps(Input.schema(), indent=4)
print(pretty_schema)

#errors = input_schema.try_set_mapping(
#    source_schema=another_input_schema,
#    function_id=uuid4(),
#    path_from=['d'],
#    path_to=['d']
#)
#print(errors)