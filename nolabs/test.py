import datetime

from pydantic.dataclasses import dataclass


@dataclass
class Test:
    value: int
    value2: datetime.datetime


t = Test('asd')