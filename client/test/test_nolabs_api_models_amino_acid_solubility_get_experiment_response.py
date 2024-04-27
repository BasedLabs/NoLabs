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

from nolabs_microservice.models.nolabs_api_models_amino_acid_solubility_get_experiment_response import NolabsApiModelsAminoAcidSolubilityGetExperimentResponse

class TestNolabsApiModelsAminoAcidSolubilityGetExperimentResponse(unittest.TestCase):
    """NolabsApiModelsAminoAcidSolubilityGetExperimentResponse unit test stubs"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def make_instance(self, include_optional) -> NolabsApiModelsAminoAcidSolubilityGetExperimentResponse:
        """Test NolabsApiModelsAminoAcidSolubilityGetExperimentResponse
            include_option is a boolean, when False only required
            params are included, when True both required and
            optional params are included """
        # uncomment below to create an instance of `NolabsApiModelsAminoAcidSolubilityGetExperimentResponse`
        """
        model = NolabsApiModelsAminoAcidSolubilityGetExperimentResponse()
        if include_optional:
            return NolabsApiModelsAminoAcidSolubilityGetExperimentResponse(
                experiment_id = None,
                experiment_name = None,
                amino_acids = None,
                properties = nolabs_microservice.models.experiment_properties_response.ExperimentPropertiesResponse(
                    amino_acid_sequence = null, 
                    fastas = null, )
            )
        else:
            return NolabsApiModelsAminoAcidSolubilityGetExperimentResponse(
                experiment_id = None,
                experiment_name = None,
                amino_acids = None,
                properties = nolabs_microservice.models.experiment_properties_response.ExperimentPropertiesResponse(
                    amino_acid_sequence = null, 
                    fastas = null, ),
        )
        """

    def testNolabsApiModelsAminoAcidSolubilityGetExperimentResponse(self):
        """Test NolabsApiModelsAminoAcidSolubilityGetExperimentResponse"""
        # inst_req_only = self.make_instance(include_optional=False)
        # inst_req_and_optional = self.make_instance(include_optional=True)

if __name__ == '__main__':
    unittest.main()