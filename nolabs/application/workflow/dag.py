import uuid
from typing import List, Type, Optional, Dict, Any
from uuid import UUID

from airflow import DAG
from airflow.decorators import task_group
from airflow.models import MappedOperator

from nolabs.application.workflow.workflow import Component
from nolabs.infrastructure.environment import Environment
from nolabs.infrastructure.settings import Settings


class DagExecutor:
    _settings: Settings

    def __init__(self):
        self._settings = Settings.load()

    async def execute(self, workflow_id: UUID, experiment_id: UUID, components: List[Component]):
        dag = generate_workflow_dag(workflow_id=workflow_id, experiment_id=experiment_id, components=components)

        if self._settings.get_environment() == Environment.LOCAL:
            dag.test()


def build_task_name(component: Component, cls: Type) -> str:
    return f'{component.name}-{cls.__name__}-{str(component.id)}'


def generate_workflow_dag(workflow_id: uuid.UUID, experiment_id: uuid.UUID, components: List[Component]):
    with DAG(
            dag_id=str(workflow_id),
            default_args={"owner": "airflow"},
            schedule_interval=None
    ) as dag:
        def component_task_group_factory(component: Component):
            group_id = f'{component.name}-{str(component.id)}'

            @task_group(group_id=group_id)
            def component_task_group():
                setup_task_name = build_task_name(component=component, cls=component.setup_operator_type)

                setup_task = component.setup_operator_type(workflow_id=workflow_id,
                                                           experiment_id=experiment_id,
                                                           component_id=component.id,
                                                           task_id=setup_task_name,
                                                           dag=dag)

                job_task: Optional[MappedOperator] = None

                if component.job_operator_type:
                    job_task_name = build_task_name(
                        component=component,
                        cls=component.job_operator_type)
                    job_task = component.job_operator_type.partial(task_id=job_task_name,
                                                                   workflow_id=workflow_id,
                                                                   experiment_id=experiment_id,
                                                                   component_id=component.id).expand(
                        job_id=setup_task.output)

                output_task_name = build_task_name(component=component, cls=component.output_operator_type)
                output_task = component.output_operator_type(component_id=component.id, experiment_id=experiment_id,
                                                             workflow_id=workflow_id, task_id=output_task_name,
                                                             dag=dag)

                if job_task:
                    job_task >> output_task
                else:
                    setup_task >> output_task

            return component_task_group()

        task_groups: Dict[uuid.UUID, Any] = {}

        for component in components:
            task_groups[component.id] = component_task_group_factory(component=component)

        # Traverse dependencies
        for component in components:
            previous_groups = [
                task_groups[prev_component_id]
                for prev_component_id in component.previous_component_ids
            ]

            if previous_groups:
                tg = task_groups[component.id]
                previous_groups >> tg  # it's ok

    return dag
