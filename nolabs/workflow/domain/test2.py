import asyncio
from typing import Optional

from pydantic import BaseModel

from nolabs.workflow.function import PythonFunction, WorkflowContext





class Input1(BaseModel):
    a: int

class SubOutput1(BaseModel):
    suboutput_test: Optional[int]

class Output1(BaseModel):
    a: int
    suboutput: SubOutput1


class Input2(BaseModel):
    b: int


class Output2(BaseModel):
    b: int


class Input3(BaseModel):
    input_param_1: int
    input_param_2: int
    input_param_3: int


class Output3(BaseModel):
    c: int


async def input_plus_one(param: Input1) -> Output1:
    return Output1(
        a=param.a + 1,
        suboutput=SubOutput1(suboutput_test=1)
    )


async def input_plus_two(param: Input2) -> Output2:
    return Output2(
        b=param.b + 2
    )


async def sum_inputs(param: Input3) -> Output3:
    return Output3(
        c=param.input_param_1 + param.input_param_2 + param.input_param_3
    )


input_plus_one_func = PythonFunction(
    pointer=input_plus_one
)
input_plus_one_func.input.set_value(value=Input1(a=1))

input_plus_two_func = PythonFunction(
    pointer=input_plus_two
)
input_plus_two_func.input.set_value(value=Input2(b=2))

sum_inputs_func = PythonFunction(
    pointer=sum_inputs
)

sum_inputs_func.set_previous([input_plus_one_func, input_plus_two_func])
sum_inputs_func.map_parameter(input_plus_one_func, ['a'], ['input_param_1'])
sum_inputs_func.map_parameter(input_plus_two_func, ['b'], ['input_param_2'])
sum_inputs_func.map_parameter(input_plus_one_func, ['suboutput','suboutput_test'], ['input_param_3'])

workflow = WorkflowContext(
    graph=[sum_inputs_func, input_plus_one_func, input_plus_two_func]
)

loop = asyncio.get_event_loop()
tasks = [
    loop.create_task(workflow.execute(terminate=True)),
]
loop.run_until_complete(asyncio.wait(tasks))
loop.close()
print('hello there')
