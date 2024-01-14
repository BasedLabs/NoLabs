# coding: utf-8

"""
    Conformations api

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 0.1.0
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


from __future__ import annotations
import json
import pprint
import re  # noqa: F401
from enum import Enum



try:
    from typing import Self
except ImportError:
    from typing_extensions import Self


class OpenMmWaterForceFields(str, Enum):
    """
    OpenMmWaterForceFields
    """

    """
    allowed enum values
    """
    AMBER14_SLASH_TIP3P_DOT_XML = 'amber14/tip3p.xml'
    AMBER14_SLASH_TIP3PFB_DOT_XML = 'amber14/tip3pfb.xml'
    AMBER14_SLASH_TIP4PEW_DOT_XML = 'amber14/tip4pew.xml'
    AMBER14_SLASH_TIP4PFB_DOT_XML = 'amber14/tip4pfb.xml'
    AMBER14_SLASH_SPCE_DOT_XML = 'amber14/spce.xml'
    AMBER14_SLASH_OPC_DOT_XML = 'amber14/opc.xml'
    AMBER14_SLASH_OPC3_DOT_XML = 'amber14/opc3.xml'
    CHARMM36_SLASH_WATER_DOT_XML = 'charmm36/water.xml'
    CHARMM36_SLASH_SPCE_DOT_XML = 'charmm36/spce.xml'
    CHARMM36_SLASH_TIP3P_MINUS_PME_MINUS_B_DOT_XML = 'charmm36/tip3p-pme-b.xml'
    CHARMM36_SLASH_TIP3P_MINUS_PME_MINUS_F_DOT_XML = 'charmm36/tip3p-pme-f.xml'
    CHARMM36_SLASH_TIP4PEW_DOT_XML = 'charmm36/tip4pew.xml'
    CHARMM36_SLASH_TIP4P2005_DOT_XML = 'charmm36/tip4p2005.xml'
    CHARMM36_SLASH_TIP5P_DOT_XML = 'charmm36/tip5p.xml'
    CHARMM36_SLASH_TIP5PEW_DOT_XML = 'charmm36/tip5pew.xml'

    @classmethod
    def from_json(cls, json_str: str) -> Self:
        """Create an instance of OpenMmWaterForceFields from a JSON string"""
        return cls(json.loads(json_str))


