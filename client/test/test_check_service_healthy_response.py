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

from nolabs_microservice.models.check_service_healthy_response import CheckServiceHealthyResponse

class TestCheckServiceHealthyResponse(unittest.TestCase):
    """CheckServiceHealthyResponse unit test stubs"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def make_instance(self, include_optional) -> CheckServiceHealthyResponse:
        """Test CheckServiceHealthyResponse
            include_option is a boolean, when False only required
            params are included, when True both required and
            optional params are included """
        # uncomment below to create an instance of `CheckServiceHealthyResponse`
        """
        model = CheckServiceHealthyResponse()
        if include_optional:
            return CheckServiceHealthyResponse(
                is_healthy = None
            )
        else:
            return CheckServiceHealthyResponse(
                is_healthy = None,
        )
        """

    def testCheckServiceHealthyResponse(self):
        """Test CheckServiceHealthyResponse"""
        # inst_req_only = self.make_instance(include_optional=False)
        # inst_req_and_optional = self.make_instance(include_optional=True)

if __name__ == '__main__':
    unittest.main()
