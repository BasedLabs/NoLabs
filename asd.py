from typing import Generic, TypeVar, reveal_type, get_args

from pydantic import BaseModel


class TestModel(BaseModel):
    number: int

TInput = TypeVar('TInput', bound=BaseModel)
TOutput = TypeVar('TOutput', bound=BaseModel)


class Test(Generic[TInput, TOutput]):
    a = TInput
    b = TOutput

    def __init__(self):
        print(self.__orig_bases__)
        print(get_args(self.__orig_bases__[0]))


class Test2(Test[int, str]):
    ...


t = Test[TestModel, TestModel]()
print(t)