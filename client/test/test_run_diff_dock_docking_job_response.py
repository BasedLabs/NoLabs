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

from nolabs_microservice.models.run_diff_dock_docking_job_response import RunDiffDockDockingJobResponse

class TestRunDiffDockDockingJobResponse(unittest.TestCase):
    """RunDiffDockDockingJobResponse unit test stubs"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def make_instance(self, include_optional) -> RunDiffDockDockingJobResponse:
        """Test RunDiffDockDockingJobResponse
            include_option is a boolean, when False only required
            params are included, when True both required and
            optional params are included """
        # uncomment below to create an instance of `RunDiffDockDockingJobResponse`
        """
        model = RunDiffDockDockingJobResponse()
        if include_optional:
            return RunDiffDockDockingJobResponse(
                predicted_pdb = None,
                predicted_ligands = None
            )
        else:
            return RunDiffDockDockingJobResponse(
                predicted_pdb = None,
                predicted_ligands = None,
        )
        """

    def testRunDiffDockDockingJobResponse(self):
        """Test RunDiffDockDockingJobResponse"""
        # inst_req_only = self.make_instance(include_optional=False)
        # inst_req_and_optional = self.make_instance(include_optional=True)

if __name__ == '__main__':
    unittest.main()
