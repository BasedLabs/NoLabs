import uuid
from collections import Counter
from typing import Optional, List
from uuid import UUID

from nolabs.exceptions import NoLabsException, ErrorCodes
from nolabs.refined.domain.models.common import Experiment
from nolabs.workflow.application.api_models import GetComponentJobIdsResponse, GetComponentJobIdsRequest, \
    GetComponentParametersResponse, GetComponentParametersRequest
from nolabs.workflow.application.models import WorkflowSchemaDbModel
from nolabs.workflow.component_factory import PythonComponentFactory
from nolabs.workflow.properties import Property
from nolabs.workflow.workflow import Workflow
from nolabs.workflow.workflow_schema import WorkflowSchemaModel, ComponentModel, PropertyModel


def map_property(property: Property) -> PropertyModel:
    return PropertyModel(
        type=property.type,
        properties={name: map_property(prop) for name, prop in (property.properties if
                                                                property.properties else {}).items()},
        required=property.required,
        description=property.description,
        enum=property.enum,
        const=property.const,
        format=property.format,
        default=property.default,
        example=property.example,
        title=property.title,
        anyOf=[
            map_property(prop)
            for prop in property.anyOf
        ],
        ref=property.ref
    )


class CreateWorkflowSchemaFeature:
    factory: PythonComponentFactory

    def __init__(self, factory: PythonComponentFactory):
        self.factory = factory

    async def handle(self, experiment_id: UUID) -> WorkflowSchemaModel:
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        workflows: WorkflowSchemaDbModel = WorkflowSchemaDbModel.objects(experiment=experiment)

        if workflows:
            return workflows[0]

        components_models: List[ComponentModel] = []

        names = Counter([comp.name for comp in self.factory.all_components()])

        for component_name in names:
            component = self.factory.create_component_instance(name=component_name)
            components_models.append(ComponentModel(
                name=component_name,
                input={name: map_property(prop) for name, prop in component.input_properties.items()},
                output={name: map_property(prop) for name, prop in component.output_properties.items()}
            ))

        workflow_schema = WorkflowSchemaModel(
            experiment_id=experiment_id,
            error=None,
            components=components_models,
            workflow_components=[]
        )

        db_model = WorkflowSchemaDbModel.create(
            id=uuid.uuid4(),
            experiment=experiment,
            value=workflow_schema
        )
        db_model.save(cascade=True)

        return workflow_schema


class GetWorkflowSchemaFeature:
    async def handle(self, experiment_id: UUID) -> Optional[WorkflowSchemaModel]:
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)

        if not db_models:
            return None

        return db_models[0].get_workflow_value()


class SetWorkflowSchemaFeature:
    factory: PythonComponentFactory

    def __init__(self, factory: PythonComponentFactory):
        self.factory = factory

    async def handle(self, workflow_schema: WorkflowSchemaModel) -> WorkflowSchemaModel:
        experiment: Experiment = Experiment.objects.with_id(workflow_schema.experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        Workflow.set_schema_errors(workflow_schema=workflow_schema, component_factory=self.factory)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)

        if not db_models:
            db_model = WorkflowSchemaDbModel.create(
                id=uuid.uuid4(),
                experiment=experiment,
                value=workflow_schema
            )
        else:
            db_model = db_models[0]
            db_model.set_workflow_value(value=workflow_schema)

        db_model.save(cascade=True)

        return workflow_schema


class StartWorkflowFeature:
    factory: PythonComponentFactory

    def __init__(self, factory: PythonComponentFactory):
        self.factory = factory

    async def handle(self, experiment_id: UUID):
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)

        if not db_models:
            raise NoLabsException(ErrorCodes.workflow_not_found)

        workflow_schema_model = db_models[0].get_workflow_value()

        workflow = Workflow(
            workflow_schema=workflow_schema_model,
            factory=self.factory
        )

        await workflow.execute()


class StopWorkflowFeature:
    factory: PythonComponentFactory

    def __init__(self, factory: PythonComponentFactory):
        self.factory = factory

    async def handle(self, experiment_id: UUID):
        experiment: Experiment = Experiment.objects.with_id(experiment_id)

        if not experiment:
            raise NoLabsException(ErrorCodes.experiment_not_found)

        db_models: List[WorkflowSchemaDbModel] = WorkflowSchemaDbModel.objects(experiment=experiment)

        if not db_models:
            raise NoLabsException(ErrorCodes.workflow_not_found)

        workflow_schema_model = db_models[0].get_workflow_value()

        workflow = Workflow(
            workflow_schema=workflow_schema_model,
            factory=self.factory
        )

        await workflow.st()


class GetComponentJobIdsFeature:
    factory: PythonComponentFactory

    def __init__(self, factory: PythonComponentFactory):
        self.factory = factory

    async def handle(self, request: GetComponentJobIdsRequest) -> GetComponentJobIdsResponse:
        if not self.factory.component_exists(name=request.name):
            raise NoLabsException(ErrorCodes.component_not_found)

        component = self.factory.create_component_instance(name=request.name, id=request.component_id)

        return GetComponentJobIdsResponse(
            job_ids=component.job_ids
        )


class GetComponentParametersFeature:
    factory: PythonComponentFactory

    def __init__(self, factory: PythonComponentFactory):
        self.factory = factory

    async def handle(self, request: GetComponentParametersRequest) -> GetComponentParametersResponse:
        if not self.factory.component_exists(name=request.name):
            raise NoLabsException(ErrorCodes.component_not_found)

        component = self.factory.create_component_instance(name=request.name, id=request.component_id)

        return GetComponentParametersResponse(
            input_dict=component.input_dict,
            output_dict=component.output_dict
        )