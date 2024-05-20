import uuid
from typing import Optional
from unittest import IsolatedAsyncioTestCase

from pydantic import BaseModel, create_model

from nolabs.workflow.component import PythonFunction
from nolabs.workflow.component import PythonComponent
from nolabs.workflow.workflow import Workflow


class TestWorkflow(IsolatedAsyncioTestCase):
    def shortDescription(self):
        return Workflow.__name__

    async def test_two_components_run(self):
        """
        Happy path
        """

        # arrange

        Input = create_model('Input', number=(int, ...))
        Output = create_model('Output', number=(int, ...))

        class ComponentNumberOne(PythonFunction[Input, Output]):
            async def execute(self):
                self.set_output_parameter(Output(number=self.input_parameter.number))

        class ComponentSimplePipe(PythonFunction[Input, Output]):
            async def execute(self):
                self.set_output_parameter(Output(number=self.input_parameter.number + 1))

        component1 = PythonComponent(
            function=ComponentNumberOne(uuid.uuid4())
        )

        component2 = PythonComponent(
            function=ComponentSimplePipe(uuid.uuid4())
        )

        component1.set_input(instance=Input(number=10))

        component2.add_previous(component1)

        component2.try_map_property(component=component1, path_from=['number'], path_to=['number'])
        workflow = Workflow(
            components=[component1, component2]
        )

        # act

        await workflow.execute()

        # assert

        self.assertEqual(component2.validate_output(), [])
        instance: Output = component2.get_output()  # type: ignore
        self.assertEqual(instance.number, 11)

    async def test_reverse_tree_components_run(self):
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

        class ComponentNumberOne(PythonFunction[Input, Output]):
            async def execute(self):
                self.set_output_parameter(Output(number=self.input_parameter.number))

        class ComponentSimplePipe(PythonFunction[Input2, Output]):
            async def execute(self):
                input_parameter = self.input_parameter
                self.set_output_parameter(Output(number=input_parameter.number1 + input_parameter.number2))

        component1 = PythonComponent(
            function=ComponentNumberOne(uuid.uuid4())
        )

        component2 = PythonComponent(
            function=ComponentNumberOne(uuid.uuid4())
        )

        component3 = PythonComponent(
            function=ComponentSimplePipe(uuid.uuid4())
        )

        component1.set_input(instance=Input(number=10))
        component2.set_input(instance=Input(number=10))
        component3.add_previous([component1, component2])
        component3.try_map_property(component=component1, path_from=['number'], path_to=['number1'])
        component3.try_map_property(component=component2, path_from=['number'], path_to=['number2'])

        workflow = Workflow(
            components=[component1, component2, component3]
        )

        # act

        await workflow.execute()

        # assert

        self.assertEqual(component3.validate_output(), [])
        instance: Output = component3.get_output()  # type: ignore
        self.assertEqual(instance.number, 20)

    async def test_tree_components_run(self):
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

        class ComponentNumberOne(PythonFunction[Input, Output]):
            async def execute(self):
                self.set_output_parameter(Output(number=self.input_parameter.number + 5))

        component1 = PythonComponent(
            function=ComponentNumberOne(uuid.uuid4())
        )

        component2 = PythonComponent(
            function=ComponentNumberOne(uuid.uuid4())
        )

        component3 = PythonComponent(
            function=ComponentNumberOne(uuid.uuid4())
        )

        component1.set_input(instance=Input(number=10))

        component2.add_previous([component1])
        component2.try_map_property(component=component1, path_from=['number'], path_to=['number'])

        component3.add_previous([component1])
        component3.try_map_property(component=component1, path_from=['number'], path_to=['number'])

        workflow = Workflow(
            components=[component1, component2, component3]
        )

        # act

        await workflow.execute()

        # assert

        self.assertEqual(component2.validate_output(), [])
        instance: Output = component2.get_output()  # type: ignore
        self.assertEqual(instance.number, 20)

        self.assertEqual(component3.validate_output(), [])
        instance: Output = component3.get_output()  # type: ignore
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

        class ComponentNumberOne(PythonFunction[Input, Output]):
            async def execute(self):
                self.set_output_parameter(
                    Output(number=self.input_parameter.number + 5)
                )

        class ComponentNumberTwo(PythonFunction[Input2, Output]):
            async def execute(self):
                self.set_output_parameter(
                    Output(
                        number=self.input_parameter.number1 + self.input_parameter.number2
                    )
                )

        component1 = PythonComponent(
            function=ComponentNumberOne(uuid.uuid4())
        )

        component2 = PythonComponent(
            function=ComponentNumberOne(uuid.uuid4())
        )

        component3 = PythonComponent(
            function=ComponentNumberOne(uuid.uuid4())
        )

        component4 = PythonComponent(
            function=ComponentNumberTwo(uuid.uuid4())
        )

        component1.set_input(instance=Input(number=10))

        component2.add_previous([component1])
        component2.try_map_property(component=component1, path_from=['number'], path_to=['number'])

        component3.add_previous([component1])
        component3.try_map_property(component=component1, path_from=['number'], path_to=['number'])

        component4.add_previous([component2, component3])
        component4.try_map_property(component=component2, path_from=['number'], path_to=['number1'])
        component4.try_map_property(component=component3, path_from=['number'], path_to=['number2'])

        workflow = Workflow(
            components=[component1, component2, component3, component4]
        )

        # act

        await workflow.execute()

        # assert

        self.assertEqual(component4.validate_output(), [])
        instance: Output = component4.get_output()  # type: ignore
        self.assertEqual(instance.number, 40)

    async def test_presevses_last_exception(self):
        """Happy path"""

        class Input(BaseModel):
            number: Optional[int] = 10

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentNumberOne(PythonFunction[Input, Output]):
            async def execute(self):
                raise Exception('Hello from the future most powerful biotech company!')  # type: ignore

        component1 = PythonComponent(
            function=ComponentNumberOne(uuid.uuid4())
        )

        workflow = Workflow(
            components=[component1]
        )

        # assert

        await workflow.execute()

        self.assertIsNotNone(component1.last_exception)
