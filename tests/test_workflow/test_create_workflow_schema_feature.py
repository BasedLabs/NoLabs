from typing import Dict, Any
from unittest import IsolatedAsyncioTestCase

from pydantic import BaseModel, create_model

from nolabs.exceptions import NoLabsException
from nolabs.workflow.application.use_cases import CreateWorkflowSchemaFeature
from nolabs.workflow.component import Component
from tests.tests_preparations import mongo_disconnect, mongo_connect, seed_experiment


class TestCreateWorkflowSchemaFeature(IsolatedAsyncioTestCase):
    def shortDescription(self) -> str:
        return CreateWorkflowSchemaFeature.__name__

    @classmethod
    def setUpClass(cls):
        mongo_connect()

    @classmethod
    def tearDownClass(cls):
        mongo_disconnect()

    async def test_cannot_register_components_with_same_ids(self):
        """
        Happy path
        """

        # arrange

        Input = create_model('Input', number=(int, ...))
        Output = create_model('Output', number=(int, ...))

        class ComponentTest(Component):
            @property
            def id(self) -> str:
                return 'hello'

            async def start(self, parameter: Input) -> Output:
                return Output(number=5)

        experiment = seed_experiment()
        sut = CreateWorkflowSchemaFeature()

        # arrange

        with self.assertRaises(NoLabsException) as context:
            workflow_schema = await sut.start(experiment_id=experiment.id, components=[
                ComponentTest(),
                ComponentTest()
            ])

    async def test_returns_components(self):
        """
        Happy path
        """

        # arrange

        Input = create_model('Input', number=(int, ...))
        Output = create_model('Output', number=(int, ...))

        class ComponentTest(Component):
            @property
            def id(self) -> str:
                return 'hello'

            async def start(self, parameter: Input) -> Output:
                return Output(number=5)

        experiment = seed_experiment()
        sut = CreateWorkflowSchemaFeature()

        # act

        workflow_schema = await sut.start(experiment_id=experiment.id, components=[
            ComponentTest()
        ])

        # arrange

        self.assertTrue(len(workflow_schema.components) > 0)

    async def test_input_output_type(self):
        """
        Happy path
        """

        class CustomModel(BaseModel):
            something: int  # type: ignore

        for _type, prop_type, prop_format, ref in [
            (int, 'integer', None, None),
            (float, 'number', None, None),
            (bool, 'boolean', None, None),
            (str, 'string', None, None),
            (bytes, 'string', 'binary', None),
            (Dict[str, Any], 'object', None, None),
            (CustomModel, None, None, '#/$defs/CustomModel')
        ]:
            # arrange

            with self.subTest(
                    msg=f'Testing that BaseModel parameter {_type} is equal to {prop_type}, {prop_format}, {ref}',
                    _type=_type,
                    prop_type=prop_type,
                    prop_format=prop_format):
                Input = create_model(
                    'Input', a=(_type, ...)
                )

                Output = create_model(
                    'Output', a=(_type, ...)
                )

                class ComponentTest(Component):
                    @property
                    def id(self) -> str:
                        return 'Test'

                    async def start(self, parameter: Input) -> Output:
                        return Output(
                            a='Something'
                        )

                sut = CreateWorkflowSchemaFeature()
                experiment = seed_experiment()

                # act

                workflow_schema = await sut.start(experiment_id=experiment.id, components=[
                    ComponentTest()
                ])

                # assert

                component = workflow_schema.components[0]

                self.assertTrue('a' in component.input)
                self.assertEqual(component.input['a'].type, prop_type)
                self.assertEqual(component.input['a'].format, prop_format)
                self.assertEqual(component.input['a'].ref, ref)

                self.assertTrue('a' in component.output)
                self.assertEqual(component.output['a'].type, prop_type)
                self.assertEqual(component.output['a'].format, prop_format)
                self.assertEqual(component.output['a'].ref, ref)

    async def test_creates_workflow(self):
        """
        Happy path
        """

        # arrange

        sut = CreateWorkflowSchemaFeature()
        experiment = seed_experiment()

        class Input(BaseModel):
            number: int  # type: ignore

        class Output(BaseModel):
            number: int  # type: ignore

        class ComponentTest(Component):

            @property
            def id(self) -> str:
                return 'Test'

            async def start(self, parameter: Input) -> Output:
                return Output(
                    number=1
                )

        # act

        workflow_schema = await sut.start(experiment_id=experiment.id, components=[
            ComponentTest()
        ])

        # assert

        self.assertIsNotNone(workflow_schema)
        self.assertIsNotNone(workflow_schema.components[0].input)
        self.assertIsNotNone(workflow_schema.components[0].output)

        properties = workflow_schema.components[0].input.items()

        self.assertTrue('number' in workflow_schema.components[0].input)
        self.assertEqual(workflow_schema.components[0].input['number'].type, 'integer')
