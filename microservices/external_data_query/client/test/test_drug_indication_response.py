# coding: utf-8

"""
    External Databases Query API

    No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)

    The version of the OpenAPI document: 0.1.0
    Generated by OpenAPI Generator (https://openapi-generator.tech)

    Do not edit the class manually.
"""  # noqa: E501


import unittest

from external_data_query_microservice.models.drug_indication_response import DrugIndicationResponse

class TestDrugIndicationResponse(unittest.TestCase):
    """DrugIndicationResponse unit test stubs"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def make_instance(self, include_optional) -> DrugIndicationResponse:
        """Test DrugIndicationResponse
            include_option is a boolean, when False only required
            params are included, when True both required and
            optional params are included """
        # uncomment below to create an instance of `DrugIndicationResponse`
        """
        model = DrugIndicationResponse()
        if include_optional:
            return DrugIndicationResponse(
                drugs = [
                    external_data_query_microservice.models.molecule.Molecule(
                        chembl_id = '', 
                        molecule_type = '', 
                        synonyms = [
                            ''
                            ], 
                        smiles = '', 
                        link = '', 
                        pref_name = '', )
                    ],
                total_count = 56,
                page = 56,
                pages = 56
            )
        else:
            return DrugIndicationResponse(
                drugs = [
                    external_data_query_microservice.models.molecule.Molecule(
                        chembl_id = '', 
                        molecule_type = '', 
                        synonyms = [
                            ''
                            ], 
                        smiles = '', 
                        link = '', 
                        pref_name = '', )
                    ],
                total_count = 56,
                page = 56,
                pages = 56,
        )
        """

    def testDrugIndicationResponse(self):
        """Test DrugIndicationResponse"""
        # inst_req_only = self.make_instance(include_optional=False)
        # inst_req_and_optional = self.make_instance(include_optional=True)

if __name__ == '__main__':
    unittest.main()