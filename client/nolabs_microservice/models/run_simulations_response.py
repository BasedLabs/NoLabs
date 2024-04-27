# coding: utf-8

"""
    NoLabs

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 1
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


from __future__ import annotations
import pprint
import re  # noqa: F401
import json


from typing import Any, ClassVar, Dict, List, Optional
from pydantic import BaseModel
from nolabs_microservice.models.pdb_content import PdbContent
try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

class RunSimulationsResponse(BaseModel):
    """
    RunSimulationsResponse
    """ # noqa: E501
    experiment_id: Optional[Any]
    experiment_name: Optional[Any]
    timeline: Optional[Any]
    pdb_content: Optional[PdbContent] = None
    __properties: ClassVar[List[str]] = ["experiment_id", "experiment_name", "timeline", "pdb_content"]

    model_config = {
        "populate_by_name": True,
        "validate_assignment": True,
        "protected_namespaces": (),
    }


    def to_str(self) -> str:
        """Returns the string representation of the model using alias"""
        return pprint.pformat(self.model_dump(by_alias=True))

    def to_json(self) -> str:
        """Returns the JSON representation of the model using alias"""
        # TODO: pydantic v2: use .model_dump_json(by_alias=True, exclude_unset=True) instead
        return json.dumps(self.to_dict())

    @classmethod
    def from_json(cls, json_str: str) -> Self:
        """Create an instance of RunSimulationsResponse from a JSON string"""
        return cls.from_dict(json.loads(json_str))

    def to_dict(self) -> Dict[str, Any]:
        """Return the dictionary representation of the model using alias.

        This has the following differences from calling pydantic's
        `self.model_dump(by_alias=True)`:

        * `None` is only added to the output dict for nullable fields that
          were set at model initialization. Other fields with value `None`
          are ignored.
        """
        _dict = self.model_dump(
            by_alias=True,
            exclude={
            },
            exclude_none=True,
        )
        # override the default output from pydantic by calling `to_dict()` of pdb_content
        if self.pdb_content:
            _dict['pdb_content'] = self.pdb_content.to_dict()
        # set to None if experiment_id (nullable) is None
        # and model_fields_set contains the field
        if self.experiment_id is None and "experiment_id" in self.model_fields_set:
            _dict['experiment_id'] = None

        # set to None if experiment_name (nullable) is None
        # and model_fields_set contains the field
        if self.experiment_name is None and "experiment_name" in self.model_fields_set:
            _dict['experiment_name'] = None

        # set to None if timeline (nullable) is None
        # and model_fields_set contains the field
        if self.timeline is None and "timeline" in self.model_fields_set:
            _dict['timeline'] = None

        return _dict

    @classmethod
    def from_dict(cls, obj: Dict) -> Self:
        """Create an instance of RunSimulationsResponse from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "experiment_id": obj.get("experiment_id"),
            "experiment_name": obj.get("experiment_name"),
            "timeline": obj.get("timeline"),
            "pdb_content": PdbContent.from_dict(obj.get("pdb_content")) if obj.get("pdb_content") is not None else None
        })
        return _obj

