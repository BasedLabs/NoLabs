from typing import Optional
from unittest import IsolatedAsyncioTestCase

from pydantic import BaseModel

from nolabs.workflow.component import Component
from nolabs.workflow.function import PythonFunction


class TestPythonFunction(IsolatedAsyncioTestCase):
    def shortDescription(self):
        return PythonFunction.__name__

    async def test_two_functions_simple_mapping(self):
        """
        Happy path
        """

        # arrange

        class Input(BaseModel):
            number: int  # type: ignore

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def start(self, parameter: Input) -> Output:
                return Output(number=parameter.number)

        class ComponentSimplePipe(Component):
            async def start(self, parameter: Input) -> Output:
                return Output(number=parameter.number)

            @property
            def name(self) -> str:
                return 'Simple pipe'

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        function2 = PythonFunction(
            component=ComponentSimplePipe()
        )

        function2.set_previous([function1])

        # act

        errors = [
            function2.try_map_property(function=function1, path_from=['number'], path_to=['number'])
        ]

        # assert

        self.assertEqual([err for err in errors if err], [])

    async def test_returns_error_if_unmapped_property_exists(self):
        # arrange
        class Input(BaseModel):
            number: int  # type: ignore

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def start(self, parameter: Input) -> Output:
                return Output(number=parameter.number + 5)

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        # act

        unmapped_properties = function1.unmapped_properties

        # assert

        self.assertEqual(len(unmapped_properties), 1)

    async def test_does_not_return_error_for_default_unmapped_property(self):
        """Happy path"""

        # arrange

        class Input(BaseModel):
            number: Optional[int] = 10

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def start(self, parameter: Input) -> Output:
                return Output(number=parameter.number + 5)  # type: ignore

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        # act

        unmapped_properties = function1.unmapped_properties

        # assert

        self.assertEqual(len(unmapped_properties), 0)

    async def test_returns_error_for_format_mismatch(self):
        # arrange
        class Input(BaseModel):
            number: Optional[int] = 10  # type: ignore

        class Output(BaseModel):
            binary: bytes  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def start(self, parameter: Input) -> Output:
                return Output(binary='hello'.encode('utf-8'))  # type: ignore

        class ComponentNumberTwo(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def start(self, parameter: Input) -> Output:
                return Output(binary='hello'.encode('utf-8'))  # type: ignore

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        function2 = PythonFunction(
            component=ComponentNumberTwo()
        )

        function2.set_previous([function1])

        # act

        error = function2.try_map_property(function1, ['binary'], ['number'])

        # assert
        self.assertIsNotNone(error)

    async def test_maps_compatible_types_int_float(self):
        """
        Happy path
        """

        # arrange

        class Input(BaseModel):
            number: Optional[int] = 10  # type: ignore

        class Output(BaseModel):
            number2: float  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def start(self, parameter: Input) -> Output:
                return Output(number2=10.0)  # type: ignore

        class ComponentNumberTwo(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def start(self, parameter: Input) -> Output:
                return Output(number2=10.0)  # type: ignore

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        function2 = PythonFunction(
            component=ComponentNumberTwo()
        )

        function2.set_previous([function1])

        # act

        error = function2.try_map_property(function1, ['number2'], ['number'])

        # assert

        self.assertIsNone(error)
