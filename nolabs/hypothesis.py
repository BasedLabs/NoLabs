import datetime
import uuid
from typing import List
from uuid import UUID

from jinja2.nodes import Test
from mongoengine import Document, ObjectIdField, UUIDField, IntField, DateTimeField, StringField, ReferenceField, \
    CASCADE, ListField, PULL
from mongoengine.base import fields
from pydantic.dataclasses import dataclass


@dataclass
class CustomInt:
    value: int

    def __int__(self) -> int:
        return self.value


class CustomIntField(fields.BaseField):

    def to_mongo(self, value):
        if isinstance(value, CustomInt):
            return value.value
        else:
            return value

    def to_python(self, value):
        if isinstance(value, int):
            return CustomInt(value)
        return value


@dataclass
class CustomId:
    value: UUID

    def __post_init__(self):
        print('Post init in id')


class CustomIdField(fields.BaseField):
    def to_mongo(self, value):
        if isinstance(value, CustomId):
            return value.value
        else:
            return value

    def to_python(self, value):
        if isinstance(value, UUID):
            return CustomId(value)
        return value.value


@dataclass
class CustomString:
    value: str

    def __str__(self):
        return self.value


class CustomString2(CustomString):
    value: str

    def __str__(self):
        return self.value


class CustomStringField(fields.BaseField):
    def to_mongo(self, value):
        if isinstance(value, CustomString):
            return value.value
        else:
            return value

    def to_python(self, value):
        if isinstance(value, str):
            return CustomString(value)
        return value


from mongoengine import connect

connect(host="mongodb://127.0.0.1:27017/my_db")


class Test2(Document):
    id: UUID = UUIDField(primary_key=True, default=uuid.uuid4)
    name: str = StringField()
    test: 'Test' = ReferenceField('Test')


class Test(Document):
    id: UUID = UUIDField(primary_key=True, default=uuid.uuid4)
    test1: CustomInt = CustomIntField()
    test2: CustomString = CustomStringField()
    tests: List[Test2] = ListField(ReferenceField(Test2))

    def __init__(self, test1: CustomInt, test2: CustomString, *args, **kwargs):
        super().__init__(test1=test1, test2=test2, *args, **kwargs)

    @property
    def iid(self) -> CustomId:
        return CustomId(self.id)


t = Test(id=uuid.uuid4(),
         test1=CustomInt(15),
         test2=CustomString('asd'))
t.save(cascade=True)
test2 = Test2(id=uuid.uuid4(),
              name='asd',
              test=t)
test2.save(cascade=True)

t = Test.objects.with_id(t.iid.value)

print(t)
