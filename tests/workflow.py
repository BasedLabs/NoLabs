import asyncio
from functools import singledispatch
from time import sleep
from typing import Optional
from unittest import IsolatedAsyncioTestCase

from pydantic import BaseModel

from nolabs.workflow.component import Component
from nolabs.workflow.function import PythonFunction
from nolabs.workflow.workflow import Workflow


class WorkflowTests(IsolatedAsyncioTestCase):
    async def test_two_functions_simple_mapping(self):
        """
        HAPPYPATH
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

            async def handle(self, parameter: Input) -> Output:
                return Output(number=parameter.number)

        class ComponentSimplePipe(Component):
            async def handle(self, parameter: Input) -> Output:
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

    async def test_two_functions_run(self):
        """
        HAPPYPATH
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

            async def handle(self, parameter: Input) -> Output:
                return Output(number=parameter.number)

        class ComponentSimplePipe(Component):
            async def handle(self, parameter: Input) -> Output:
                return Output(number=parameter.number + 1)

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

        function1.set_input(instance=Input(number=10))
        function2.try_map_property(function=function1, path_from=['number'], path_to=['number'])
        workflow = Workflow(
            functions=[function1, function2]
        )

        # act

        await workflow.execute()

        # assert

        self.assertEqual(function2.validate_output(), [])
        instance: Output = function2.get_output()  # type: ignore
        self.assertEqual(instance.number, 11)

    async def test_reverse_tree_functions_run(self):
        """
        HAPPYPATH
        O
            > O
        O

        O is node
        > shows data flow vector
        """

        class Input(BaseModel):
            number: int  # type: ignore

        class Output(BaseModel):
            number: int  # type: ignore

        class Input2(BaseModel):
            number1: int  # type: ignore
            number2: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                return Output(number=parameter.number)

        class ComponentSimplePipe(Component):
            async def handle(self, parameter: Input2) -> Output:
                return Output(number=parameter.number1 + parameter.number2)

            @property
            def name(self) -> str:
                return 'Simple pipe'

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        function2 = PythonFunction(
            component=ComponentNumberOne()
        )

        function3 = PythonFunction(
            component=ComponentSimplePipe()
        )

        function1.set_input(instance=Input(number=10))
        function2.set_input(instance=Input(number=10))
        function3.set_previous([function1, function2])
        function3.try_map_property(function=function1, path_from=['number'], path_to=['number1'])
        function3.try_map_property(function=function2, path_from=['number'], path_to=['number2'])

        workflow = Workflow(
            functions=[function1, function2, function3]
        )

        # act

        await workflow.execute()

        # assert

        self.assertEqual(function3.validate_output(), [])
        instance: Output = function3.get_output()  # type: ignore
        self.assertEqual(instance.number, 20)

    async def test_tree_functions_run(self):
        """
        HAPPYPATH
            O
        O >
            O
        O is node
        > shows data flow vector
        """

        class Input(BaseModel):
            number: int  # type: ignore

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                return Output(number=parameter.number + 5)


        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        function2 = PythonFunction(
            component=ComponentNumberOne()
        )

        function3 = PythonFunction(
            component=ComponentNumberOne()
        )

        function1.set_input(instance=Input(number=10))

        function2.set_previous([function1])
        function2.try_map_property(function=function1, path_from=['number'], path_to=['number'])

        function3.set_previous([function1])
        function3.try_map_property(function=function1, path_from=['number'], path_to=['number'])

        workflow = Workflow(
            functions=[function1, function2, function3]
        )

        # act

        await workflow.execute()

        # assert

        self.assertEqual(function2.validate_output(), [])
        instance: Output = function2.get_output()  # type: ignore
        self.assertEqual(instance.number, 20)

        self.assertEqual(function3.validate_output(), [])
        instance: Output = function3.get_output()  # type: ignore
        self.assertEqual(instance.number, 20)

    async def test_diamond_run(self):
        """
        HAPPYPATH
            O
        O >   > O
            O
        O is node
        > shows data flow vector
        """

        class Input(BaseModel):
            number: int  # type: ignore

        class Output(BaseModel):
            number: int  # type: ignore

        class Input2(BaseModel):
            number1: int  # type: ignore
            number2: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                return Output(number=parameter.number + 5)

        class ComponentNumberTwo(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input2) -> Output:
                return Output(number=parameter.number1 + parameter.number2)


        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        function2 = PythonFunction(
            component=ComponentNumberOne()
        )

        function3 = PythonFunction(
            component=ComponentNumberOne()
        )

        function4 = PythonFunction(
            component=ComponentNumberTwo()
        )

        function1.set_input(instance=Input(number=10))

        function2.set_previous([function1])
        function2.try_map_property(function=function1, path_from=['number'], path_to=['number'])

        function3.set_previous([function1])
        function3.try_map_property(function=function1, path_from=['number'], path_to=['number'])

        function4.set_previous([function2, function3])
        function4.try_map_property(function=function2, path_from=['number'], path_to=['number1'])
        function4.try_map_property(function=function3, path_from=['number'], path_to=['number2'])

        workflow = Workflow(
            functions=[function1, function2, function3, function4]
        )

        # act

        await workflow.execute()

        # assert

        self.assertEqual(function4.validate_output(), [])
        instance: Output = function4.get_output()  # type: ignore
        self.assertEqual(instance.number, 40)

    async def test_returns_error_if_unmapped_property_exists(self):

        class Input(BaseModel):
            number: int  # type: ignore

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                return Output(number=parameter.number + 5)

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        # assert

        unmapped_properties = function1.unmapped_properties
        self.assertEqual(len(unmapped_properties), 1)

    async def test_does_not_return_error_for_default_unmapped_property(self):
        """HAPPYPATH"""

        class Input(BaseModel):
            number: Optional[int] = 10  # type: ignore

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                return Output(number=parameter.number + 5)  # type: ignore

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        # assert

        unmapped_properties = function1.unmapped_properties
        self.assertEqual(len(unmapped_properties), 0)

    async def test_returns_error_for_format_mismatch(self):
        class Input(BaseModel):
            number: Optional[int] = 10  # type: ignore

        class Output(BaseModel):
            binary: bytes  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                return Output(binary='hello'.encode('utf-8'))  # type: ignore

        class ComponentNumberTwo(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                return Output(binary='hello'.encode('utf-8'))  # type: ignore

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        function2 = PythonFunction(
            component=ComponentNumberTwo()
        )

        function2.set_previous([function1])

        # assert

        error = function2.try_map_property(function1, ['binary'], ['number'])
        self.assertIsNotNone(error)

    async def test_maps_compatible_types_int_float(self):
        """
        HAPPYPATH
        """
        class Input(BaseModel):
            number: Optional[int] = 10  # type: ignore

        class Output(BaseModel):
            number2: float  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                return Output(number2=10.0)  # type: ignore

        class ComponentNumberTwo(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                return Output(number2=10.0)  # type: ignore

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        function2 = PythonFunction(
            component=ComponentNumberTwo()
        )

        function2.set_previous([function1])

        # assert

        error = function2.try_map_property(function1, ['number2'], ['number'])
        self.assertIsNone(error)

    async def test_presevses_last_exception(self):
        """HAPPYPATH"""

        class Input(BaseModel):
            number: Optional[int] = 10  # type: ignore

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                raise Exception('Hello from the future most powerful biotech company!')  # type: ignore

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        workflow = Workflow(
            functions=[function1]
        )

        # assert

        await workflow.execute()

        self.assertIsNotNone(function1.last_exception)

    async def test_terminates_long_running(self):
        """HAPPYPATH"""

        class Input(BaseModel):
            number: Optional[int] = 10  # type: ignore

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def handle(self, parameter: Input) -> Output:
                sleep(1000)
                raise Exception('Hello from the future most powerful biotech company!')  # type: ignore

        function1 = PythonFunction(
            component=ComponentNumberOne()
        )

        workflow = Workflow(
            functions=[function1]
        )

        loop = asyncio.get_event_loop()

        async def wrapper(delay, coro):
            await asyncio.sleep(delay)
            return await coro

        # assert

        loop.create_task(workflow.execute())
        await wrapper(2.0, function1.terminate())

