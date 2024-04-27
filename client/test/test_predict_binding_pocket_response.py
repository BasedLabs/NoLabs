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

from nolabs_microservice.models.predict_binding_pocket_response import PredictBindingPocketResponse

class TestPredictBindingPocketResponse(unittest.TestCase):
    """PredictBindingPocketResponse unit test stubs"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def make_instance(self, include_optional) -> PredictBindingPocketResponse:
        """Test PredictBindingPocketResponse
            include_option is a boolean, when False only required
            params are included, when True both required and
            optional params are included """
        # uncomment below to create an instance of `PredictBindingPocketResponse`
        """
        model = PredictBindingPocketResponse()
        if include_optional:
            return PredictBindingPocketResponse(
                pocket_ids = None
            )
        else:
            return PredictBindingPocketResponse(
                pocket_ids = None,
        )
        """

    def testPredictBindingPocketResponse(self):
        """Test PredictBindingPocketResponse"""
        # inst_req_only = self.make_instance(include_optional=False)
        # inst_req_and_optional = self.make_instance(include_optional=True)

if __name__ == '__main__':
    unittest.main()