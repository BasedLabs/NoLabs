import uuid
from typing import List
from unittest import IsolatedAsyncioTestCase

from pydantic import BaseModel, create_model
from workflow import (DefaultWorkflowComponentModelValue, MappingModel,
                      WorkflowComponentModel)

from nolabs.application.use_cases.workflow.use_cases import (
    CreateWorkflowSchemaFeature, UpdateWorkflowSchemaFeature)
from tests.test_workflow.mixins import WorkflowTestsMixin
from tests.tests_preparations import mongo_connect


class TestUpdateWorkflowSchemaFeature(IsolatedAsyncioTestCase, WorkflowTestsMixin):
    def setUp(self):
        mongo_connect()

    async def test_component_not_exists_error(self):
        # arrange

        Input = create_model("Input", **{"a": (int, ...)})
        Output = create_model("Output", **{"b": (float, ...)})
        C = self.seed_empty_component(Input, Output)

        create_schema_feature = CreateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="Wf1",
                component_id=uuid.uuid4(),
                mappings=[],
                error=None,
                defaults=[],
            )
        )

        # act

        schema = await set_schema_feature.handle(workflow_schema=schema)

        # assert

        self.assertFalse(schema.valid)
        self.assertTrue(schema.workflow_components[0].error)

    async def test_component_mapping_one_level_primitive(self):
        # arrange

        Input = create_model("Input", **{"a": (int, ...)})
        Output = create_model("Output", **{"b": (int, ...)})
        C = self.seed_empty_component(Input, Output)

        create_schema_feature = CreateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        comp_id_1 = uuid.uuid4()
        comp_id_2 = uuid.uuid4()
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C", component_id=comp_id_1, mappings=[], error=None, defaults=[]
            )
        )
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C",
                component_id=comp_id_2,
                mappings=[
                    MappingModel(
                        source_path=["b"],
                        target_path=["a"],
                        source_component_id=comp_id_1,
                    )
                ],
                error=None,
                defaults=[],
            )
        )

        # act

        schema = await set_schema_feature.handle(workflow_schema=schema)

        # assert

        self.assertTrue(schema.valid)
        self.assertFalse(schema.workflow_components[1].mappings[0].error)

    async def test_component_mapping_one_level_list(self):
        # arrange

        Input = create_model("Input", **{"a": (List[int], ...)})
        Output = create_model("Output", **{"b": (List[float], ...)})
        C = self.seed_empty_component(Input, Output)

        create_schema_feature = CreateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        comp_id_1 = uuid.uuid4()
        comp_id_2 = uuid.uuid4()
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C", component_id=comp_id_1, mappings=[], error=None, defaults=[]
            )
        )
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C",
                component_id=comp_id_2,
                mappings=[
                    MappingModel(
                        source_path=["b"],
                        target_path=["a"],
                        source_component_id=comp_id_1,
                    )
                ],
                error=None,
                defaults=[],
            )
        )

        # act

        schema = await set_schema_feature.handle(workflow_schema=schema)

        # assert

        self.assertTrue(schema.valid)
        self.assertFalse(schema.workflow_components[1].mappings[0].error)

    async def test_component_mapping_two_levels(self):
        # arrange

        class InnerInput(BaseModel):
            b: int

        class Input(BaseModel):
            a: InnerInput

        Output = create_model("Output", **{"b": (int, ...)})
        C = self.seed_empty_component(Input, Output)

        create_schema_feature = CreateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        comp_id_1 = uuid.uuid4()
        comp_id_2 = uuid.uuid4()
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C", component_id=comp_id_1, mappings=[], error=None, defaults=[]
            )
        )
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C",
                component_id=comp_id_2,
                mappings=[
                    MappingModel(
                        source_path=["b"],
                        target_path=["a", "b"],
                        source_component_id=comp_id_1,
                    )
                ],
                error=None,
                defaults=[],
            )
        )

        # act

        schema = await set_schema_feature.handle(workflow_schema=schema)

        # assert

        self.assertTrue(schema.valid)
        self.assertFalse(schema.workflow_components[1].mappings[0].error)

    async def test_component_mapping_incorrect_path(self):
        # arrange

        Input = create_model("Input", **{"a": (int, ...)})

        Output = create_model("Output", **{"b": (int, ...)})
        C = self.seed_empty_component(Input, Output)

        create_schema_feature = CreateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        comp_id_1 = uuid.uuid4()
        comp_id_2 = uuid.uuid4()
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C", component_id=comp_id_1, mappings=[], error=None, defaults=[]
            )
        )
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C",
                component_id=comp_id_2,
                mappings=[
                    MappingModel(
                        source_path=["ccc"],
                        target_path=["ccc"],
                        source_component_id=comp_id_1,
                    )
                ],
                error=None,
                defaults=[],
            )
        )

        # act

        schema = await set_schema_feature.handle(workflow_schema=schema)

        # assert

        self.assertFalse(schema.valid)
        self.assertTrue(schema.workflow_components[1].mappings[0].error)

    async def test_component_property_default(self):
        # arrange

        Input = create_model("Input", **{"a": (int, ...)})

        Output = create_model("Output", **{"b": (int, ...)})
        C = self.seed_empty_component(Input, Output)

        create_schema_feature = CreateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        comp_id_1 = uuid.uuid4()
        comp_id_2 = uuid.uuid4()
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C", component_id=comp_id_1, mappings=[], error=None, defaults=[]
            )
        )
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C",
                component_id=comp_id_2,
                mappings=[
                    MappingModel(
                        source_path=["b"],
                        target_path=["a"],
                        source_component_id=comp_id_1,
                    )
                ],
                error=None,
                defaults=[
                    DefaultWorkflowComponentModelValue(target_path=["a"], value=15)
                ],
            )
        )

        # act

        schema = await set_schema_feature.handle(workflow_schema=schema)

        # assert

        self.assertTrue(schema.valid)

    async def test_component_property_default_incompatible_types(self):
        # arrange

        Input = create_model("Input", **{"a": (bytes, ...)})

        Output = create_model("Output", **{"b": (int, ...)})
        C = self.seed_empty_component(Input, Output)

        create_schema_feature = CreateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        comp_id_1 = uuid.uuid4()
        comp_id_2 = uuid.uuid4()
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C", component_id=comp_id_1, mappings=[], error=None, defaults=[]
            )
        )
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C",
                component_id=comp_id_2,
                mappings=[],
                error=None,
                defaults=[
                    DefaultWorkflowComponentModelValue(target_path=["a"], value=15)
                ],
            )
        )

        # act

        schema = await set_schema_feature.handle(workflow_schema=schema)

        # assert

        self.assertFalse(schema.valid)
        self.assertTrue(schema.workflow_components[1].defaults[0].error)

    async def test_component_property_default_compatible_types(self):
        # arrange

        Input = create_model("Input", **{"a": (int, ...)})

        Output = create_model("Output", **{"b": (float, ...)})
        C = self.seed_empty_component(Input, Output)

        create_schema_feature = CreateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        comp_id_1 = uuid.uuid4()
        comp_id_2 = uuid.uuid4()
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C", component_id=comp_id_1, mappings=[], error=None, defaults=[]
            )
        )
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C",
                component_id=comp_id_2,
                mappings=[],
                error=None,
                defaults=[
                    DefaultWorkflowComponentModelValue(target_path=["a"], value=15)
                ],
            )
        )

        # act

        schema = await set_schema_feature.handle(workflow_schema=schema)

        # assert

        self.assertTrue(schema.valid)
        self.assertFalse(schema.workflow_components[1].defaults[0].error)

    async def test_component_property_default_inner_types(self):
        # arrange

        Input2 = create_model("Input2", **{"a": (int, ...)})

        Input = create_model("Input", **{"c": (Input2, ...)})

        Output = create_model("Output", **{"b": (float, ...)})
        C = self.seed_empty_component(Input, Output)

        create_schema_feature = CreateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        set_schema_feature = UpdateWorkflowSchemaFeature(
            available_components={C.name: C}
        )

        experiment = self.seed_experiment()

        schema = await create_schema_feature.handle(experiment_id=experiment.iid.value)
        comp_id_1 = uuid.uuid4()
        comp_id_2 = uuid.uuid4()
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C", component_id=comp_id_1, mappings=[], error=None, defaults=[]
            )
        )
        schema.workflow_components.append(
            WorkflowComponentModel(
                name="C",
                component_id=comp_id_2,
                mappings=[],
                error=None,
                defaults=[
                    DefaultWorkflowComponentModelValue(target_path=["c", "a"], value=15)
                ],
            )
        )

        # act

        schema = await set_schema_feature.handle(workflow_schema=schema)

        # assert

        self.assertTrue(schema.valid)
        self.assertFalse(schema.workflow_components[1].defaults[0].error)
