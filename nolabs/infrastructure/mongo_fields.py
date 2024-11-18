from mongoengine import fields

from nolabs.seedwork.domain.value_objects import ValueObjectFloat, ValueObjectString


class ValueObjectStringField(fields.BaseField):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if "factory" not in kwargs:
            raise ValueError

        self.factory = kwargs["factory"]

    def to_mongo(self, value):
        if isinstance(value, ValueObjectString):
            return value.value
        else:
            return value

    def to_python(self, value):
        if isinstance(value, str):
            return self.factory(value)
        return value


class ValueObjectFloatField(fields.BaseField):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if "factory" not in kwargs:
            raise ValueError

        self.factory = kwargs["factory"]

    def to_mongo(self, value):
        if isinstance(value, ValueObjectFloat):
            return value.value
        else:
            return value

    def to_python(self, value):
        if isinstance(value, float):
            return self.factory(value)
        return value
