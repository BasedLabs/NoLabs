from typing import Type

from biobuddy_microservice import DefaultApi
from pydantic import BaseModel

from nolabs.workflow.component import PythonComponent


class InputA(BaseModel):
    number: int


class OutputA(BaseModel):
    number: int


class InputPlusOneComponent(PythonComponent[InputA, OutputA]):
    @property
    def _input_parameter_type(self) -> Type[InputA]:
        return InputA

    @property
    def _output_parameter_type(self) -> Type[OutputA]:
        return OutputA

    api: DefaultApi
    name = 'Input plus one'

    def _execute(self):
        inp = self.input
        output = OutputA(
            number=inp.number + 1
        )
        self.output = output

    async def _restore_parameters(self):
        pass


class InputPlusTwoComponent(PythonComponent[InputA, OutputA]):
    @property
    def _input_parameter_type(self) -> Type[InputA]:
        return InputA

    @property
    def _output_parameter_type(self) -> Type[OutputA]:
        return OutputA

    api: DefaultApi
    name = 'Input plus two'

    def _execute(self):
        inp = self.input
        output = OutputA(
            number=inp.number + 2
        )
        self.output = output

    async def _restore_parameters(self):
        pass