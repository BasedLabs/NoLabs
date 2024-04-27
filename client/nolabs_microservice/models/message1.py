# coding: utf-8

"""
    NoLabs

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 1
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


from __future__ import annotations
from inspect import getfullargspec
import json
import pprint
import re  # noqa: F401

from typing import Optional
from pydantic import BaseModel, Field, StrictStr, ValidationError, field_validator
from nolabs_microservice.models.function_call import FunctionCall
from nolabs_microservice.models.regular_message import RegularMessage
from typing import Union, Any, List, TYPE_CHECKING, Optional, Dict
from typing_extensions import Literal
from pydantic import StrictStr, Field
try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

MESSAGE1_ANY_OF_SCHEMAS = ["FunctionCall", "RegularMessage"]

class Message1(BaseModel):
    """
    Message1
    """

    # data type: RegularMessage
    anyof_schema_1_validator: Optional[RegularMessage] = None
    # data type: FunctionCall
    anyof_schema_2_validator: Optional[FunctionCall] = None
    if TYPE_CHECKING:
        actual_instance: Optional[Union[FunctionCall, RegularMessage]] = None
    else:
        actual_instance: Any = None
    any_of_schemas: List[str] = Literal[MESSAGE1_ANY_OF_SCHEMAS]

    model_config = {
        "validate_assignment": True,
        "protected_namespaces": (),
    }

    def __init__(self, *args, **kwargs) -> None:
        if args:
            if len(args) > 1:
                raise ValueError("If a position argument is used, only 1 is allowed to set `actual_instance`")
            if kwargs:
                raise ValueError("If a position argument is used, keyword arguments cannot be used.")
            super().__init__(actual_instance=args[0])
        else:
            super().__init__(**kwargs)

    @field_validator('actual_instance')
    def actual_instance_must_validate_anyof(cls, v):
        instance = Message1.model_construct()
        error_messages = []
        # validate data type: RegularMessage
        if not isinstance(v, RegularMessage):
            error_messages.append(f"Error! Input type `{type(v)}` is not `RegularMessage`")
        else:
            return v

        # validate data type: FunctionCall
        if not isinstance(v, FunctionCall):
            error_messages.append(f"Error! Input type `{type(v)}` is not `FunctionCall`")
        else:
            return v

        if error_messages:
            # no match
            raise ValueError("No match found when setting the actual_instance in Message1 with anyOf schemas: FunctionCall, RegularMessage. Details: " + ", ".join(error_messages))
        else:
            return v

    @classmethod
    def from_dict(cls, obj: dict) -> Self:
        return cls.from_json(json.dumps(obj))

    @classmethod
    def from_json(cls, json_str: str) -> Self:
        """Returns the object represented by the json string"""
        instance = cls.model_construct()
        error_messages = []
        # anyof_schema_1_validator: Optional[RegularMessage] = None
        try:
            instance.actual_instance = RegularMessage.from_json(json_str)
            return instance
        except (ValidationError, ValueError) as e:
             error_messages.append(str(e))
        # anyof_schema_2_validator: Optional[FunctionCall] = None
        try:
            instance.actual_instance = FunctionCall.from_json(json_str)
            return instance
        except (ValidationError, ValueError) as e:
             error_messages.append(str(e))

        if error_messages:
            # no match
            raise ValueError("No match found when deserializing the JSON string into Message1 with anyOf schemas: FunctionCall, RegularMessage. Details: " + ", ".join(error_messages))
        else:
            return instance

    def to_json(self) -> str:
        """Returns the JSON representation of the actual instance"""
        if self.actual_instance is None:
            return "null"

        to_json = getattr(self.actual_instance, "to_json", None)
        if callable(to_json):
            return self.actual_instance.to_json()
        else:
            return json.dumps(self.actual_instance)

    def to_dict(self) -> Dict:
        """Returns the dict representation of the actual instance"""
        if self.actual_instance is None:
            return "null"

        to_json = getattr(self.actual_instance, "to_json", None)
        if callable(to_json):
            return self.actual_instance.to_dict()
        else:
            return json.dumps(self.actual_instance)

    def to_str(self) -> str:
        """Returns the string representation of the actual instance"""
        return pprint.pformat(self.model_dump())

