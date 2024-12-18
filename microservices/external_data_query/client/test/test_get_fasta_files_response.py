# coding: utf-8

"""
    RCSB PDB Query API

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 0.1.0
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


import unittest

from external_data_query_microservice.models.get_fasta_files_response import (
    GetFastaFilesResponse,
)


class TestGetFastaFilesResponse(unittest.TestCase):
    """GetFastaFilesResponse unit test stubs"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def make_instance(self, include_optional) -> GetFastaFilesResponse:
        """Test GetFastaFilesResponse
        include_option is a boolean, when False only required
        params are included, when True both required and
        optional params are included"""
        # uncomment below to create an instance of `GetFastaFilesResponse`
        """
        model = GetFastaFilesResponse()
        if include_optional:
            return GetFastaFilesResponse(
                fasta_contents = [
                    external_data_query_microservice.models.fetched_protein.FetchedProtein(
                        fasta_contents = '', 
                        link = '', )
                    ]
            )
        else:
            return GetFastaFilesResponse(
                fasta_contents = [
                    external_data_query_microservice.models.fetched_protein.FetchedProtein(
                        fasta_contents = '', 
                        link = '', )
                    ],
        )
        """

    def testGetFastaFilesResponse(self):
        """Test GetFastaFilesResponse"""
        # inst_req_only = self.make_instance(include_optional=False)
        # inst_req_and_optional = self.make_instance(include_optional=True)


if __name__ == "__main__":
    unittest.main()
