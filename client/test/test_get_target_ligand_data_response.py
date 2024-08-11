# coding: utf-8

"""
    NoLabs

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 1
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


import unittest
import datetime

from nolabs_microservice.models.get_target_ligand_data_response import GetTargetLigandDataResponse

class TestGetTargetLigandDataResponse(unittest.TestCase):
    """GetTargetLigandDataResponse unit test stubs"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def make_instance(self, include_optional) -> GetTargetLigandDataResponse:
        """Test GetTargetLigandDataResponse
            include_option is a boolean, when False only required
            params are included, when True both required and
            optional params are included """
        # uncomment below to create an instance of `GetTargetLigandDataResponse`
        """
        model = GetTargetLigandDataResponse()
        if include_optional:
            return GetTargetLigandDataResponse(
                ligand_id = None,
                ligand_name = None,
                ligand_sdf = None,
                ligand_smiles = None
            )
        else:
            return GetTargetLigandDataResponse(
                ligand_id = None,
                ligand_name = None,
                ligand_sdf = None,
                ligand_smiles = None,
        )
        """

    def testGetTargetLigandDataResponse(self):
        """Test GetTargetLigandDataResponse"""
        # inst_req_only = self.make_instance(include_optional=False)
        # inst_req_and_optional = self.make_instance(include_optional=True)

if __name__ == '__main__':
    unittest.main()