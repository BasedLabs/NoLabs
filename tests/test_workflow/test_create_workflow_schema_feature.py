from typing import List
from unittest import IsolatedAsyncioTestCase

from pydantic import create_model, BaseModel

from nolabs.application.use_cases.workflow.use_cases import CreateWorkflowSchemaFeature
from tests.test_workflow.mixins import WorkflowTestsMixin
from tests.tests_preparations import mongo_connect


class TestCreateWorkflowSchemaFeature(IsolatedAsyncioTestCase, WorkflowTestsMixin):
    def setUp(self):
        mongo_connect()

    async def test_component_exists(self):
        # arrange

        Input = create_model('Input', **{
            'a': (int,...)
        })
        Output = create_model('Output', **{
            'b': (float,...)
        })

        C = self.seed_empty_component(Input, Output)

        feature = CreateWorkflowSchemaFeature({
            C.name: C
        })

        experiment = self.seed_experiment()

        # act

        schema = await feature.handle(experiment_id=experiment.iid.value)

        # assert

        self.assertEqual(len(schema.components), 1)

    async def test_primitive_type_io(self):
        # arrange

        Input = create_model('Input', **{
            'a': (int,...)
        })
        Output = create_model('Output', **{
            'b': (float,...)
        })

        C = self.seed_empty_component(Input, Output)

        feature = CreateWorkflowSchemaFeature({
            C.name: C
        })

        experiment = self.seed_experiment()

        # act

        schema = await feature.handle(experiment_id=experiment.iid.value)
        component = schema.components[0]

        # assert

        self.assertEqual(component.input['a'].type, 'integer')
        self.assertEqual(component.output['b'].type, 'number')

    async def test_complex_type(self):
        # arrange

        class InnerInput(BaseModel):
            i: int = 10

        class Input(BaseModel):
            a: InnerInput

        Output = create_model('Output', **{
            'b': (float,...)
        })

        C = self.seed_empty_component(Input, Output)

        feature = CreateWorkflowSchemaFeature({
            C.name: C
        })

        experiment = self.seed_experiment()

        # act

        schema = await feature.handle(experiment_id=experiment.iid.value)
        component = schema.components[0]

        # assert

        self.assertEqual(component.input['a'].properties['i'].type, 'integer')

    async def test_list_property(self):
        # arrange

        class Input(BaseModel):
            a: List[int]

        Output = create_model('Output', **{
            'b': (float,...)
        })

        C = self.seed_empty_component(Input, Output)

        feature = CreateWorkflowSchemaFeature({
            C.name: C
        })

        experiment = self.seed_experiment()

        # act

        schema = await feature.handle(experiment_id=experiment.iid.value)
        component = schema.components[0]

        # assert

        self.assertEqual(component.input['a'].items.type, 'integer')

    async def test_default_value(self):
        # arrange

        class Input(BaseModel):
            a: int = 10

        Output = create_model('Output', **{
            'b': (float,...)
        })

        C = self.seed_empty_component(Input, Output)

        feature = CreateWorkflowSchemaFeature({
            C.name: C
        })

        experiment = self.seed_experiment()

        # act

        schema = await feature.handle(experiment_id=experiment.iid.value)
        component = schema.components[0]

        # assert

        self.assertEqual(component.input['a'].default, 10)
