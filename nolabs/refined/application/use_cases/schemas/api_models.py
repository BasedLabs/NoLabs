from enum import Enum
from typing import List

from pydantic.dataclasses import dataclass


@dataclass
class RequestParameter:
    name: str
    type: str


@dataclass
class Response:
    name: str
    type: str


class EndpointMethod(str, Enum):
    post = 'POST'
    get = 'GET'
    delete = 'DELETE'
    patch = 'PATCH'


@dataclass
class Endpoint:
    path: str
    method: EndpointMethod

    request_parameters: List[RequestParameter]
    response: Response


@dataclass
class ApiSchema:
    tag: str

    endpoints: List[Endpoint]

