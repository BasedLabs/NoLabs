import uuid
from typing import Optional, Type
from unittest import IsolatedAsyncioTestCase

from pydantic import create_model

from nolabs.workflow.component import PythonComponent, TOutput, TInput


class TestPythonComponents(IsolatedAsyncioTestCase):
    def shortDescription(self):
        return PythonComponent.__name__

    async def test_two_components_simple_mapping(self):
        """
        Happy path
        """

        # arrange

        Input = create_model('Input', number=(int, ...))
        Output = create_model('Output', number=(int, ...))

        class PythonNumberOne(PythonComponent[Input, Output]):

            async def restore_output(self):
                pass

            @property
            def _input_parameter_type(self) -> Type[TInput]:
                pass

            @property
            def _output_parameter_type(self) -> Type[TOutput]:
                pass

            async def _execute(self):
                self.output = Output(number=10)

        class PythonSimplePipe(PythonComponent[Input, Output]):
            async def _execute(self):
                self.output = Output(
                    number=self.input.number + 1
                )

        component1 = PythonNumberOne(uuid.uuid4())

        component2 = PythonSimplePipe(uuid.uuid4())

        component2.add_previous(component1)

        # act

        errors = [
            component2.try_map_property(component=component1, path_from=['number'], path_to=['number'])
        ]

        # assert

        self.assertEqual([err for err in errors if err], [])

    async def test_returns_error_if_unmapped_property_exists(self):
        # arrange

        Input = create_model('Input', number=(int, ...))
        Output = create_model('Output', number=(int, ...))

        class PythonNumberOne(PythonComponent[Input, Output]):
            async def _execute(self):
                self.output = {
                    'number': 10
                }

        component1 = PythonNumberOne(uuid.uuid4())

        # act

        unmapped_properties = component1.unmapped_properties

        # assert

        self.assertEqual(len(unmapped_properties), 1)

    async def test_does_not_return_error_for_default_unmapped_property(self):
        """Happy path"""

        # arrange

        Input = create_model('Input', number=(Optional[int], 10))
        Output = create_model('Output', number=(int, ...))

        class PythonNumberOne(PythonComponent[Input, Output]):
            async def restore_output(self):
                pass

            @property
            def _input_parameter_type(self) -> Type[Input]:
                return Input

            @property
            def _output_parameter_type(self) -> Type[Output]:
                return Output

            async def _execute(self):
                self.output = Output(
                    number=10
                )

        component1 = PythonNumberOne(uuid.uuid4())

        # act

        unmapped_properties = component1.unmapped_properties

        # assert

        self.assertEqual(len(unmapped_properties), 0)

    async def test_returns_error_for_format_mismatch(self):
        # arrange
        Input = create_model('Input', number=(int, ...))
        Output = create_model('Output', binary=(bytes, ...))

        class PythonNumberOne(PythonComponent[Input, Output]):
            async def _execute(self):
                self.output = Output(
                    binary='Hello there'.encode()
                )

        class PythonNumberTwo(PythonComponent[Input, Output]):
            async def execute(self):
                self.output = Output(binary='hello'.encode('utf-8'))  # type: ignore

        component1 = PythonNumberOne(uuid.uuid4())

        component2 = PythonNumberTwo((uuid.uuid4()))

        component2.add_previous(component1)

        # act

        error = component2.try_map_property(component1, ['binary'], ['number'])

        # assert
        self.assertIsNotNone(error)

    async def test_maps_compatible_types_int_float(self):
        """
        Happy path
        """

        # arrange

        Input = create_model('Input', number=(Optional[int], 10))
        Output = create_model('Output', number2=(float, ...))

        class PythonNumberOne(PythonComponent[Input, Output]):
            async def execute(self):
                self.output = Output(number2=10.0)

        class PythonNumberTwo(PythonComponent[Input, Output]):
            async def execute(self):
                self.output = Output(number2=10.0)  # type: ignore

        component1 = PythonNumberOne(uuid.uuid4())
        component2 = PythonNumberTwo(uuid.uuid4())

        component2.add_previous(component1)

        # act

        error = component2.try_map_property(component1, ['number2'], ['number'])

        # assert

        self.assertIsNone(error)
