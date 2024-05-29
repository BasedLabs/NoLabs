from mongoengine import Document, ReferenceField, IntField

from mongoengine import connect

connect(host="mongodb://127.0.0.1:27017/my_db")

class A(Document):
    number1: int = IntField()

    meta = {
        'allow_inheritance': True
    }


class B(A):
    number2: int = IntField()


class C(Document):
    a: A = ReferenceField(A)


a = B(
        number1=20,
        number2=15
    )
a.save()
c = C(
    a=a
)
c.save()

c = C.objects.with_id(c.id)
print(c)


