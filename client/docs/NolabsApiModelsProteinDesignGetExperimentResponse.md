# NolabsApiModelsProteinDesignGetExperimentResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **object** |  | 
**experiment_name** | **object** |  | 
**pdb_files** | **object** |  | 
**properties** | [**NolabsApiModelsProteinDesignExperimentPropertiesResponse**](NolabsApiModelsProteinDesignExperimentPropertiesResponse.md) |  | 

## Example

```python
from nolabs_microservice.models.nolabs_api_models_protein_design_get_experiment_response import NolabsApiModelsProteinDesignGetExperimentResponse

# TODO update the JSON string below
json = "{}"
# create an instance of NolabsApiModelsProteinDesignGetExperimentResponse from a JSON string
nolabs_api_models_protein_design_get_experiment_response_instance = NolabsApiModelsProteinDesignGetExperimentResponse.from_json(json)
# print the JSON string representation of the object
print NolabsApiModelsProteinDesignGetExperimentResponse.to_json()

# convert the object into a dict
nolabs_api_models_protein_design_get_experiment_response_dict = nolabs_api_models_protein_design_get_experiment_response_instance.to_dict()
# create an instance of NolabsApiModelsProteinDesignGetExperimentResponse from a dict
nolabs_api_models_protein_design_get_experiment_response_form_dict = nolabs_api_models_protein_design_get_experiment_response.from_dict(nolabs_api_models_protein_design_get_experiment_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


