import datetime
import uuid
from typing import List
from uuid import UUID

from mongoengine import Document, ObjectIdField, UUIDField, IntField, DateTimeField, StringField, ReferenceField, \
    CASCADE, ListField, PULL
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


class CustomIdField(UUIDField):

    def to_python(self, value):
        if isinstance(value, CustomId):
            return value.value
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
    id: CustomId = CustomIdField(db_field='_id', primary_key=True)
    test1: CustomInt = CustomIntField()
    test2: CustomString2 = CustomStringField()

    @classmethod
    def pre_save(cls, sender, document, **kwargs):
        pass

class Test2(Document):
    id: CustomId = CustomIdField(db_field='_id', primary_key=True)
    test1: CustomInt = CustomIntField()
    test2: CustomString2 = CustomStringField()
    reference = ListField(ReferenceField(Test, reverse_delete_rule=PULL))


objects: List[Test] = list(Test.objects.all())
t = Test(id=CustomId(uuid.uuid4()),
         test1=CustomInt(10),
         test2=CustomString2('asd'))
t3 = Test(id=CustomId(uuid.uuid4()),
         test1=CustomInt(15),
         test2=CustomString2('asd2222'))
t2 = Test2(
    id=CustomId(uuid.uuid4()),
    test1=CustomInt(10),
    test2=CustomString2('asd'),
    reference=[t, t3]
)

# t = Test(id=uuid.uuid4(),
#         value=10,
#         value2='asd')
t.save()
t2.save()
t2.save(cascade=True)

t = Test.objects.get(id=t.id)
t.delete()

t2 = Test2.objects.get(id=t2.id)

print(t)
