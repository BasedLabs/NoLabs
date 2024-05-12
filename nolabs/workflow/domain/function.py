from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, List, Callable
from uuid import UUID, uuid4


@dataclass
class InputParameter:
    name: str
    type: str
    description: Optional[str]


@dataclass
class OutputParameter:
    type: str
    description: Optional[str]


class Function(ABC):
    id: UUID
    name: str
    signature: str

    inputs: List[InputParameter]
    output: OutputParameter

    parents: List['Function']
    children: List['Function']

    def __init__(self,
                 id: UUID,
                 name: str,
                 signature: str,
                 inputs: List[InputParameter],
                 output: OutputParameter,
                 parents: List['Function'],
                 children: List['Function']):
        self.id = id
        self.name = name
        self.signature = signature
        self.inputs = inputs
        self.output = output
        self.parents = parents
        self.children = children

    @abstractmethod
    def execute(self):
        pass


@dataclass
class PythonClass:
    init_dependencies: List[InputParameter]


class PythonFunction(Function):
    pointer: Callable
    _inputs_set: bool = False

    def __init__(self, pointer: Callable, id: UUID, name: str, signature: str, inputs: List[InputParameter],
                 output: OutputParameter, parents: List['Function'], children: List['Function']):

        super().__init__(id, name, signature, inputs, output, parents, children)
        self.pointer = pointer

    def _find_input_parameter(self, name: str) -> Optional[InputParameter]:
        for input_parameter in self.inputs:
            if input_parameter.name == name:
                return input_parameter

        return None

    def set_inputs(self, *args, **kwargs):
        self._inputs_set = True

        for input_name, input_value in kwargs:
            for input_parameter in self.inputs:


class WorkflowContext:
    def __init__(self, functions: List[Function]):
