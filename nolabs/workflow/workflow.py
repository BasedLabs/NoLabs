from dataclasses import dataclass
from typing import List

from nolabs.workflow.component import Component
from nolabs.workflow.models import ComponentDbModel, WorkflowSchemaDbModel
from nolabs.workflow.workflow_schema import WorkflowSchemaModel, JobValidationError


@dataclass
class WorkflowValidationError:
    msg: str


class Workflow:
    async def execute(self, workflow_schema: WorkflowSchemaModel, components: List[Component]):
        workflow_db_model: WorkflowSchemaDbModel = WorkflowSchemaDbModel.objects.with_id(workflow_schema.workflow_id)

        async def execute(component: Component):
            component_db_model: ComponentDbModel = ComponentDbModel.objects.with_id(component.id)

            for previous in component.previous:
                await execute(previous)

                validation_errors = previous.validate_output()

                if validation_errors:
                    workflow_schema.get_wf_component(previous.id).error = ', '.join(
                        [ve.msg for ve in validation_errors])
                    workflow_schema.valid = False
                    workflow_db_model.set_workflow_value(workflow_schema)
                    workflow_db_model.save()

                return

            input_changed = component.set_input_from_previous()

            if not input_changed:
                return

            component_db_model.input_parameter_dict = component.input_parameter_dict
            component_db_model.save()

            validation_errors = component.validate_input()

            if validation_errors:
                workflow_schema.get_wf_component(component.id).error = ', '.join(
                    [ve.msg for ve in validation_errors])
                workflow_db_model.set_workflow_value(workflow_schema)
                workflow_db_model.save()
                return

            await component.setup_jobs()
            component_db_model.jobs = component.jobs
            component_db_model.save()

            jobs_validation_errors = await component.prevalidate_jobs()

            if jobs_validation_errors:
                workflow_schema.get_wf_component(component.id).error = [JobValidationError(
                    job_id=ve.job_id,
                    msg=ve.msg
                ) for ve in jobs_validation_errors]
            else:
                workflow_schema.get_wf_component(component.id).jobs_errors = []
                workflow_db_model.set_workflow_value(workflow_schema)
                workflow_db_model.save()

            await component.execute()
            component_db_model.output_parameter_dict = component.output_parameter_dict
            component_db_model.save()

        for component in components:
            await execute(component=component)
