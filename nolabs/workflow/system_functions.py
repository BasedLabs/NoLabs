from typing import Any, List

from nolabs.workflow.function import PythonFunction


async def take_property(obj: Any, prop_name: str) -> Any:
    return getattr(obj, prop_name)


async def concat_lists(obj1: List[Any], obj2: List[Any]) -> List[Any]:
    if obj1 is None and obj2 is None:
        return []

    if obj1 is None:
        return obj2 if isinstance(obj2, list) else []

    if obj2 is None:
        return obj1 if isinstance(obj1, list) else []

    if not isinstance(obj1, list):
        obj1 = [obj1]

    if not isinstance(obj2, list):
        obj2 = [obj2]

    if not obj1:
        return obj2

    if not obj2:
        return obj1

    if obj1[0] is not obj2[0]:
        raise ValueError('Types do not match')

    return obj1 + obj2


def take_property_workflow_function_factory() -> PythonFunction:
    return PythonFunction(
        pointer=take_property,
    )


def concat_lists_workflow_function_factory() -> PythonFunction:
    return PythonFunction(
        pointer=concat_lists
    )
