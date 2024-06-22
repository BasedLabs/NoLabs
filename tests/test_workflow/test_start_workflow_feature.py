import uuid
from typing import Type, List
from unittest import IsolatedAsyncioTestCase

from pydantic import create_model

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.application.use_cases.workflow.use_cases import CreateWorkflowSchemaFeature, StartWorkflowFeature, \
    UpdateWorkflowSchemaFeature
from nolabs.application.workflow.component import Component, TOutput, TInput, JobValidationError
from nolabs.application.workflow.models import ComponentDbModel
from nolabs.application.workflow.workflow_schema import WorkflowComponentModel, MappingModel, DefaultWorkflowComponentModelValue
from tests.test_workflow.mixins import WorkflowTestsMixin
from tests.tests_preparations import mongo_connect


class TestUpdateWorkflowSchemaFeature(IsolatedAsyncioTestCase, WorkflowTestsMixin):
    def setUp(self):
        mongo_connect()

    async def test_linear_execute(self):
        """
        X - X
        """

        # arrange

        Input = create_model('Input', **{
            'a': (int, ...)
        })
        Output = create_model('Output', **{
            'b': (int, ...)
        })

        class C(Component[Input, Output]):
            name = 'C'

            async def setup_jobs(self):
                pass

            async def jobs_setup_errors(self) -> List[JobValidationError]:
                return []

            @property
            def _input_parameter_type(self) -> Type[TInput]:
                return Input

            @property
            def _output_parameter_type(self) -> Type[TOutput]:
                return Output

            async def execute(self):
                self.output = Output(
                    b=self.input.a
                )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={
                C.name: C
            }
        )

        create_schema_feature = CreateWorkflowSchemaFeature(available_components={
            C.name: C
        })

        start_workflow_feature = StartWorkflowFeature(available_components={
            C.name: C
        })

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        id1 = uuid.uuid4()
        id2 = uuid.uuid4()
        schema.workflow_components = [
            WorkflowComponentModel(
                name='C',
                component_id=id1,
                mappings=[],
                error=None,
                defaults=[
                    DefaultWorkflowComponentModelValue(
                        target_path=['a'],
                        value=15
                    )
                ]
            ),
            WorkflowComponentModel(
                name='C',
                component_id=id2,
                mappings=[MappingModel(
                    source_path=['b'],
                    target_path=['a'],
                    source_component_id=id1
                )],
                error=None,
                defaults=[]
            )
        ]
        await set_schema_feature.handle(workflow_schema=schema)

        # act

        await start_workflow_feature.handle(workflow_id=schema.workflow_id)

        # assert

        component_db_model: ComponentDbModel = ComponentDbModel.objects.with_id(id2)
        output = component_db_model.output_parameter_dict
        self.assertEquals(output['b'], 15)

    async def test_diamond_test_execute(self):
        # arrange

        LinearInput = create_model('LinearInput', **{
            'a': (int, ...),
        })

        DuplexInput = create_model('DuplexInput', **{
            'x1': (int, ...),
            'x2': (int, ...)
        })

        Output = create_model('Output', **{
            'b': (int, ...)
        })

        class Linear(Component[LinearInput, Output]):
            name = 'Linear'

            async def setup_jobs(self):
                pass

            async def jobs_setup_errors(self) -> List[JobValidationError]:
                return []

            @property
            def _input_parameter_type(self) -> Type[TInput]:
                return LinearInput

            @property
            def _output_parameter_type(self) -> Type[TOutput]:
                return Output

            async def execute(self):
                self.output = Output(
                    b=self.input.a + 1
                )

        class Duplex(Component[DuplexInput, Output]):
            name = 'Duplex'

            async def setup_jobs(self):
                pass

            async def jobs_setup_errors(self) -> List[JobValidationError]:
                return []

            @property
            def _input_parameter_type(self) -> Type[TInput]:
                return DuplexInput

            @property
            def _output_parameter_type(self) -> Type[TOutput]:
                return Output

            async def execute(self):
                self.output = Output(
                    b=self.input.x1 + self.input.x2
                )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={
                Linear.name: Linear,
                Duplex.name: Duplex
            }
        )

        create_schema_feature = CreateWorkflowSchemaFeature(available_components={
            Linear.name: Linear,
            Duplex.name: Duplex
        })

        start_workflow_feature = StartWorkflowFeature(available_components={
            Linear.name: Linear,
            Duplex.name: Duplex
        })

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        id1 = uuid.uuid4()
        id2 = uuid.uuid4()
        id3 = uuid.uuid4()
        id4 = uuid.uuid4()
        schema.workflow_components = [
            WorkflowComponentModel(
                name='Linear',
                component_id=id1,
                mappings=[],
                error=None,
                defaults=[
                    DefaultWorkflowComponentModelValue(
                        target_path=['a'],
                        value=1
                    )
                ]
            ),
            WorkflowComponentModel(
                name='Linear',
                component_id=id2,
                mappings=[MappingModel(
                    source_path=['b'],
                    target_path=['a'],
                    source_component_id=id1
                )],
                error=None
            ),
            WorkflowComponentModel(
                name='Linear',
                component_id=id3,
                mappings=[MappingModel(
                    source_path=['b'],
                    target_path=['a'],
                    source_component_id=id1
                )],
                error=None
            ),
            WorkflowComponentModel(
                name='Duplex',
                component_id=id4,
                mappings=[MappingModel(
                    source_path=['b'],
                    target_path=['x1'],
                    source_component_id=id2
                ),
                    MappingModel(
                        source_path=['b'],
                        target_path=['x2'],
                        source_component_id=id3
                    )
                ],
                error=None,
                defaults=[]
            )
        ]
        schema = await set_schema_feature.handle(workflow_schema=schema)

        # act

        await start_workflow_feature.handle(workflow_id=schema.workflow_id)

        # assert

        component_db_model: ComponentDbModel = ComponentDbModel.objects.with_id(id4)
        output = component_db_model.output_parameter_dict
        self.assertEquals(output['b'], 6)

    async def test_component_jobs_exception(self):
        """
        X - X
        """

        # arrange

        Input = create_model('Input', **{
            'a': (int, ...)
        })
        Output = create_model('Output', **{
            'b': (int, ...)
        })

        class C(Component[Input, Output]):
            name = 'C'

            async def setup_jobs(self):
                pass

            async def jobs_setup_errors(self) -> List[JobValidationError]:
                return [
                    JobValidationError(job_id=uuid.uuid4(), msg='Issue 1'),
                    JobValidationError(job_id=uuid.uuid4(), msg='Issue 2')
                ]

            @property
            def _input_parameter_type(self) -> Type[TInput]:
                return Input

            @property
            def _output_parameter_type(self) -> Type[TOutput]:
                return Output

            async def execute(self):
                self.output = Output(
                    b=self.input.a
                )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={
                C.name: C
            }
        )

        create_schema_feature = CreateWorkflowSchemaFeature(available_components={
            C.name: C
        })

        start_workflow_feature = StartWorkflowFeature(available_components={
            C.name: C
        })

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        id1 = uuid.uuid4()
        schema.workflow_components = [
            WorkflowComponentModel(
                name='C',
                component_id=id1,
                mappings=[],
                error=None,
                defaults=[
                    DefaultWorkflowComponentModelValue(
                        target_path=['a'],
                        value=15
                    )
                ]
            )
        ]
        await set_schema_feature.handle(workflow_schema=schema)

        # act

        await start_workflow_feature.handle(workflow_id=schema.workflow_id)

        # assert

        component: ComponentDbModel = ComponentDbModel.objects.with_id(id1)

        self.assertTrue(component.jobs_errors)

    async def test_component_exception(self):
        """
        X - X
        """

        # arrange

        Input = create_model('Input', **{
            'a': (int, ...)
        })
        Output = create_model('Output', **{
            'b': (int, ...)
        })

        class C(Component[Input, Output]):
            name = 'C'

            async def setup_jobs(self):
                pass

            async def jobs_setup_errors(self) -> List[JobValidationError]:
                return []

            @property
            def _input_parameter_type(self) -> Type[TInput]:
                return Input

            @property
            def _output_parameter_type(self) -> Type[TOutput]:
                return Output

            async def execute(self):
                raise NoLabsException(ErrorCodes.experiment_not_found)

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={
                C.name: C
            }
        )

        create_schema_feature = CreateWorkflowSchemaFeature(available_components={
            C.name: C
        })

        start_workflow_feature = StartWorkflowFeature(available_components={
            C.name: C
        })

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        id1 = uuid.uuid4()
        schema.workflow_components = [
            WorkflowComponentModel(
                name='C',
                component_id=id1,
                mappings=[],
                error=None,
                defaults=[
                    DefaultWorkflowComponentModelValue(
                        target_path=['a'],
                        value=15
                    )
                ]
            )
        ]
        await set_schema_feature.handle(workflow_schema=schema)

        # act

        await start_workflow_feature.handle(workflow_id=schema.workflow_id)

        # assert

        component: ComponentDbModel = ComponentDbModel.objects.with_id(id1)

        self.assertTrue(component.last_exceptions)

