from typing import Optional
from unittest import IsolatedAsyncioTestCase

from pydantic import BaseModel

from nolabs.workflow.component import Component
from nolabs.workflow.function import PythonFunction
from nolabs.workflow.workflow import Workflow


class TestWorkflow(IsolatedAsyncioTestCase):
    def shortDescription(self):
        return Workflow.__name__

    async def test_two_functions_run(self):
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
        Happy path
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

            async def start(self, parameter: Input) -> Output:
                return Output(number=parameter.number)

        class ComponentSimplePipe(Component):
            async def start(self, parameter: Input2) -> Output:
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
        Happy path
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

            async def start(self, parameter: Input) -> Output:
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
        Happy path
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

            async def start(self, parameter: Input) -> Output:
                return Output(number=parameter.number + 5)

        class ComponentNumberTwo(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def start(self, parameter: Input2) -> Output:
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

    async def test_presevses_last_exception(self):
        """Happy path"""

        class Input(BaseModel):
            number: Optional[int] = 10

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentNumberOne(Component):
            @property
            def name(self) -> str:
                return 'Simple pipe'

            async def start(self, parameter: Input) -> Output:
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
