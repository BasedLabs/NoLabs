from dataclasses import dataclass
from typing import List

from nolabs.exceptions import NoLabsException
from nolabs.workflow.component import Component
from nolabs.workflow.models import ComponentDbModel, WorkflowSchemaDbModel, JobErrorDbModel, InputPropertyErrorDbModel
from nolabs.workflow.workflow_schema import WorkflowSchemaModel


@dataclass
class WorkflowValidationError:
    msg: str


class Workflow:
    async def execute_single(self, component: Component):
        component_db_model: ComponentDbModel = ComponentDbModel.objects.with_id(component.id)

        for previous in component.previous:
            validation_errors = previous.validate_output()

            if validation_errors:
                component_db_model.input_property_errors = [InputPropertyErrorDbModel(loc=err.loc, msg=err.msg) for err in validation_errors]
                component_db_model.save()
                return

        input_changed = component.set_input_from_previous()

        if input_changed:
            component_db_model.input_parameter_dict = component.input_parameter_dict
            component_db_model.save()

            validation_errors = component.validate_input()

            if validation_errors:
                component_db_model.input_property_errors = [
                    InputPropertyErrorDbModel.create(
                        loc=e.loc,
                        msg=e.msg
                    ) for e in validation_errors
                ]
                component_db_model.save()
                return

            await component.setup_jobs()
            component_db_model.jobs = component.jobs
            component_db_model.save()

            jobs_validation_errors = await component.prevalidate_jobs()

            if jobs_validation_errors:
                component_db_model.jobs_errors = [JobErrorDbModel(
                    job_id=ve.job_id,
                    msg=ve.msg
                ) for ve in jobs_validation_errors]
            else:
                component_db_model.jobs_errors = []

            component_db_model.save()

        if component.validate_output():
            try:
                await component.execute()
                component_db_model.output_parameter_dict = component.output_parameter_dict
                component_db_model.save()
            except Exception as e:
                if isinstance(e, NoLabsException):
                    component_db_model.last_exceptions = e.messages
                else:
                    component_db_model.last_exceptions = [str(e)]

                component_db_model.save()

    async def execute(self, components: List[Component]):
        async def execute(component: Component):
            component_db_model: ComponentDbModel = ComponentDbModel.objects.with_id(component.id)

            for previous in component.previous:
                await execute(previous)

                validation_errors = previous.validate_output()

                if validation_errors:
                    component_db_model.input_property_errors = [InputPropertyErrorDbModel(loc=err.loc, msg=err.msg) for err in validation_errors]
                    component_db_model.save()
                    return

            input_changed = component.set_input_from_previous()

            if input_changed:
                component_db_model.input_parameter_dict = component.input_parameter_dict
                component_db_model.save()

                validation_errors = component.validate_input()

                if validation_errors:
                    component_db_model.input_property_errors = [
                        InputPropertyErrorDbModel.create(
                            loc=e.loc,
                            msg=e.msg
                        ) for e in validation_errors
                    ]
                    component_db_model.save()
                    return

                await component.setup_jobs()
                component_db_model.jobs = component.jobs
                component_db_model.save()

                jobs_validation_errors = await component.prevalidate_jobs()

                if jobs_validation_errors:
                    component_db_model.jobs_errors = [JobErrorDbModel(
                        job_id=ve.job_id,
                        msg=ve.msg
                    ) for ve in jobs_validation_errors]
                else:
                    component_db_model.jobs_errors = []

                component_db_model.save()

            if component.validate_output():
                try:
                    await component.execute()
                    component_db_model.output_parameter_dict = component.output_parameter_dict
                    component_db_model.save()
                except Exception as e:
                    if isinstance(e, NoLabsException):
                        component_db_model.last_exceptions = e.messages
                    else:
                        component_db_model.last_exceptions = [str(e)]

                    component_db_model.save()

        for component in components:
            await execute(component=component)
