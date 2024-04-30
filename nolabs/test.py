import datetime
import uuid
from typing import List
from uuid import UUID

from mongoengine import Document, ObjectIdField, UUIDField, IntField, DateTimeField, StringField, ReferenceField, \
    CASCADE, ListField, PULL
from mongoengine.base import fields
from pydantic.dataclasses import dataclass


@dataclass
class CustomInt:
    value: int

    def __int__(self) -> int:
        return self.value


class CustomIntField(IntField):

    def to_python(self, value):
        if isinstance(value, CustomInt):
            return value.value
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
        return value

@dataclass
class CustomString:
    value: str

    def __str__(self):
        return self.value


class CustomString2(CustomString):
    value: str

    def __str__(self):
        return self.value


class CustomStringField(StringField):

    def to_python(self, value):
        if isinstance(value, CustomString):
            return value.value
        return value


from mongoengine import connect

connect(host="mongodb://127.0.0.1:27017/my_db")


class Test(Document):
    id: CustomId = CustomIdField(db_field="_id", primary_key=True)
    test1: CustomInt = CustomIntField()
    test2: CustomString2 = CustomStringField()

    def __init__(self, id: CustomId, test1: CustomInt, test2: CustomString2, *args, **kwargs):
        super().__init__(id=id, test1=test1, test2=test2, *args, **kwargs)

        if int(test1) != 10:
            raise ValueError


t = Test(id=CustomId(uuid.uuid4()),
         test1=CustomInt(10),
         test2=CustomString2('asd'))
t.save()

Test.objects.filter(id=t.id.value).delete()

t2: Test = Test.objects.with_id(t.id.value)
print(t)
