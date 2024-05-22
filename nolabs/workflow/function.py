from abc import abstractmethod
from typing import TypeVar, Generic, Any, Dict, Union, Type, get_args, List
from uuid import UUID

from pydantic import BaseModel

TInput = TypeVar('TInput', bound=BaseModel)
TOutput = TypeVar('TOutput', bound=BaseModel)


class PythonFunctionState(BaseModel):
    running: bool = False


class PythonFunction(Generic[TInput, TOutput]):
    id: UUID

    _input_parameter_dict: Dict[str, Any]
    _output_parameter_dict: Dict[str, Any]

    input_parameter_type: Type[TInput]
    output_parameter_type: Type[TOutput]

    job_ids: List[UUID]

    def __init__(self, id: UUID):
        self.id = id

        self.input_parameter_type, self.output_parameter_type = get_args(self.__orig_bases__[0])

        self._input_parameter_dict = {}
        self._output_parameter_dict = {}

    async def state(self) -> PythonFunctionState:
        return PythonFunctionState(
            running=False
        )

    async def restore(self):
        ...

    async def stop(self):
        ...

    @abstractmethod
    async def execute(self):
        ...

    def set_input_parameter(self, input_parameter: Union[TInput, Dict[str, Any]]):
        if isinstance(input_parameter, BaseModel):
            self._input_parameter_dict = input_parameter.dict()
            return

        self._input_parameter_dict = input_parameter

    def set_output_parameter(self, output_parameter: Union[TOutput, Dict[str, Any]]):
        if isinstance(output_parameter, BaseModel):
            self._output_parameter_dict = output_parameter.dict()
            return

        self._output_parameter_dict = output_parameter

    @property
    def input_parameter_dict(self) -> Dict[str, Any]:
        return self._input_parameter_dict

    @property
    def output_parameter_dict(self) -> Dict[str, Any]:
        return self._output_parameter_dict

    @property
    def input_parameter(self) -> TInput:
        return self.input_parameter_type(**self._input_parameter_dict)

    @property
    def output_parameter(self) -> TOutput:
        return self.output_parameter_type(**self._output_parameter_dict)
