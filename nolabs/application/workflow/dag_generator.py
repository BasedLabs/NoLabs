import uuid
from typing import List, Type, Optional, Dict, Any

from airflow import DAG
from airflow.decorators import task_group

from nolabs.application.workflow.component import Component, JobOperator


def build_task_name(component: Component, cls: Type) -> str:
    return f'{component.name}-{cls.__name__}-{str(component.id)}'


def generate_workflow_dag(workflow_id: uuid.UUID, components: List[Component]):
    dag = DAG(
        dag_id=str(workflow_id),
        default_args={"owner": "airflow"},
        schedule_interval=None
    )

    def component_task_group_factory(component: Component):
        group_id = f'{component.name}-{str(component.id)}'

        @task_group(group_id=group_id)
        def component_task_group():
            setup_task_name = build_task_name(component=component, cls=component.setup_operator_type)

            setup_task = component.setup_operator_type(component_id=component.id, task_id=setup_task_name, dag=dag)

            job_task: Optional[JobOperator] = None

            if component.job_operator_type:
                job_task_name = build_task_name(component=component, cls=component.job_operator_type)
                component.job_operator_type.partial(task_id=job_task_name, component_id=component.id).expand(
                    job_id=setup_task.output)
                job_task = component.job_operator_type(component_id=component.id, task_id=job_task_name, dag=dag)

            output_task_name = build_task_name(component=component, cls=component.output_operator_type)
            output_task = component.output_operator_type(output_task_name, dag=dag)

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

        tg = task_groups[component.id]

        previous_groups >> tg  # it's ok

    return dag
