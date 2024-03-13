# coding: utf-8

"""
    ChemBL Query API

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 0.1.0
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


from __future__ import annotations
import pprint
import re  # noqa: F401
import json

from pydantic import BaseModel, StrictInt, StrictStr
from typing import Any, ClassVar, Dict, List, Optional
from chembl_query_microservice.models.job_id import JobId
from typing import Optional, Set
from typing_extensions import Self

class ChEMBLMoleculeRequest(BaseModel):
    """
    ChEMBLMoleculeRequest
    """ # noqa: E501
    filters: Optional[Dict[str, Any]] = None
    order_by: Optional[StrictStr] = None
    limit: Optional[StrictInt] = 20
    job_id: Optional[JobId] = None
    __properties: ClassVar[List[str]] = ["filters", "order_by", "limit", "job_id"]

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
    def from_json(cls, json_str: str) -> Optional[Self]:
        """Create an instance of ChEMBLMoleculeRequest from a JSON string"""
        return cls.from_dict(json.loads(json_str))

    def to_dict(self) -> Dict[str, Any]:
        """Return the dictionary representation of the model using alias.

        This has the following differences from calling pydantic's
        `self.model_dump(by_alias=True)`:

        * `None` is only added to the output dict for nullable fields that
          were set at model initialization. Other fields with value `None`
          are ignored.
        """
        excluded_fields: Set[str] = set([
        ])

        _dict = self.model_dump(
            by_alias=True,
            exclude=excluded_fields,
            exclude_none=True,
        )
        # override the default output from pydantic by calling `to_dict()` of job_id
        if self.job_id:
            _dict['job_id'] = self.job_id.to_dict()
        return _dict

    @classmethod
    def from_dict(cls, obj: Optional[Dict[str, Any]]) -> Optional[Self]:
        """Create an instance of ChEMBLMoleculeRequest from a dict"""
        if obj is None:
            return None

        if not isinstance(obj, dict):
            return cls.model_validate(obj)

        _obj = cls.model_validate({
            "filters": obj.get("filters"),
            "order_by": obj.get("order_by"),
            "limit": obj.get("limit") if obj.get("limit") is not None else 20,
            "job_id": JobId.from_dict(obj["job_id"]) if obj.get("job_id") is not None else None
        })
        return _obj


