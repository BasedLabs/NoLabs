# coding: utf-8

"""
    Reinvent4 API

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 0.1
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


import unittest

from reinvent_microservice.models.smiles_response import SmilesResponse

class TestSmilesResponse(unittest.TestCase):
    """SmilesResponse unit test stubs"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def make_instance(self, include_optional) -> SmilesResponse:
        """Test SmilesResponse
            include_option is a boolean, when False only required
            params are included, when True both required and
            optional params are included """
        # uncomment below to create an instance of `SmilesResponse`
        """
        model = SmilesResponse()
        if include_optional:
            return SmilesResponse(
                smiles = [
                    reinvent_microservice.models.smiles.Smiles(
                        smiles = '', 
                        drug_likeness = 1.337, 
                        score = 1.337, )
                    ]
            )
        else:
            return SmilesResponse(
                smiles = [
                    reinvent_microservice.models.smiles.Smiles(
                        smiles = '', 
                        drug_likeness = 1.337, 
                        score = 1.337, )
                    ],
        )
        """

    def testSmilesResponse(self):
        """Test SmilesResponse"""
        # inst_req_only = self.make_instance(include_optional=False)
        # inst_req_and_optional = self.make_instance(include_optional=True)

if __name__ == '__main__':
    unittest.main()