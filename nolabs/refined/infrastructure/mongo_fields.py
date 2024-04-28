from mongoengine import StringField, UUIDField, FloatField, BinaryField, FileField

from nolabs.seedwork.domain.value_objects import ValueObjectString, ValueObjectUUID, ValueObjectFloat, ValueObjectBinary


class ValueObjectStringField(StringField):
    def to_python(self, value):
        if isinstance(value, ValueObjectString):
            return value.value
        return value


class ValueObjectUUIDField(UUIDField):
    def to_python(self, value):
        if isinstance(value, ValueObjectUUID):
            return value.value
        return value


class ValueObjectFloatField(FloatField):
    def to_python(self, value):
        if isinstance(value, ValueObjectFloat):
            return value.value
        return value