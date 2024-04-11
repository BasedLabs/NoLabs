# NolabsApiModelsAminoAcidFoldingGetExperimentResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **object** |  | 
**experiment_name** | **object** |  | 
**amino_acids** | **object** |  | 
**properties** | [**NolabsApiModelsAminoAcidFoldingExperimentPropertiesResponse**](NolabsApiModelsAminoAcidFoldingExperimentPropertiesResponse.md) |  | 

## Example

```python
from nolabs_microservice.models.nolabs_api_models_amino_acid_folding_get_experiment_response import NolabsApiModelsAminoAcidFoldingGetExperimentResponse

# TODO update the JSON string below
json = "{}"
# create an instance of NolabsApiModelsAminoAcidFoldingGetExperimentResponse from a JSON string
nolabs_api_models_amino_acid_folding_get_experiment_response_instance = NolabsApiModelsAminoAcidFoldingGetExperimentResponse.from_json(json)
# print the JSON string representation of the object
print NolabsApiModelsAminoAcidFoldingGetExperimentResponse.to_json()

# convert the object into a dict
nolabs_api_models_amino_acid_folding_get_experiment_response_dict = nolabs_api_models_amino_acid_folding_get_experiment_response_instance.to_dict()
# create an instance of NolabsApiModelsAminoAcidFoldingGetExperimentResponse from a dict
nolabs_api_models_amino_acid_folding_get_experiment_response_form_dict = nolabs_api_models_amino_acid_folding_get_experiment_response.from_dict(nolabs_api_models_amino_acid_folding_get_experiment_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


