# coding: utf-8

"""
    Reinvent4 API

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 0.1
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


import unittest

from reinvent_microservice.models.response_smiles_jobs_job_id_smiles_get import ResponseSmilesJobsJobIdSmilesGet

class TestResponseSmilesJobsJobIdSmilesGet(unittest.TestCase):
    """ResponseSmilesJobsJobIdSmilesGet unit test stubs"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def make_instance(self, include_optional) -> ResponseSmilesJobsJobIdSmilesGet:
        """Test ResponseSmilesJobsJobIdSmilesGet
            include_option is a boolean, when False only required
            params are included, when True both required and
            optional params are included """
        # uncomment below to create an instance of `ResponseSmilesJobsJobIdSmilesGet`
        """
        model = ResponseSmilesJobsJobIdSmilesGet()
        if include_optional:
            return ResponseSmilesJobsJobIdSmilesGet(
                smiles = None
            )
        else:
            return ResponseSmilesJobsJobIdSmilesGet(
                smiles = None,
        )
        """

    def testResponseSmilesJobsJobIdSmilesGet(self):
        """Test ResponseSmilesJobsJobIdSmilesGet"""
        # inst_req_only = self.make_instance(include_optional=False)
        # inst_req_and_optional = self.make_instance(include_optional=True)

if __name__ == '__main__':
    unittest.main()