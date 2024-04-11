# NolabsApiModelsAminoAcidFoldingExperimentPropertiesResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**amino_acid_sequence** | [**AminoAcidSequence**](AminoAcidSequence.md) |  | 
**fastas** | **object** |  | 

## Example

```python
from nolabs_microservice.models.nolabs_api_models_amino_acid_folding_experiment_properties_response import NolabsApiModelsAminoAcidFoldingExperimentPropertiesResponse

# TODO update the JSON string below
json = "{}"
# create an instance of NolabsApiModelsAminoAcidFoldingExperimentPropertiesResponse from a JSON string
nolabs_api_models_amino_acid_folding_experiment_properties_response_instance = NolabsApiModelsAminoAcidFoldingExperimentPropertiesResponse.from_json(json)
# print the JSON string representation of the object
print NolabsApiModelsAminoAcidFoldingExperimentPropertiesResponse.to_json()

# convert the object into a dict
nolabs_api_models_amino_acid_folding_experiment_properties_response_dict = nolabs_api_models_amino_acid_folding_experiment_properties_response_instance.to_dict()
# create an instance of NolabsApiModelsAminoAcidFoldingExperimentPropertiesResponse from a dict
nolabs_api_models_amino_acid_folding_experiment_properties_response_form_dict = nolabs_api_models_amino_acid_folding_experiment_properties_response.from_dict(nolabs_api_models_amino_acid_folding_experiment_properties_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


