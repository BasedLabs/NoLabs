# NolabsApiModelsSmallMoleculesDesignGetExperimentResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **object** |  | 
**experiment_name** | **object** |  | 
**created_at** | **object** |  | 
**status** | [**GetExperimentStatusResponse**](GetExperimentStatusResponse.md) |  | 
**properties** | [**NolabsApiModelsSmallMoleculesDesignExperimentPropertiesResponse**](NolabsApiModelsSmallMoleculesDesignExperimentPropertiesResponse.md) |  | 

## Example

```python
from nolabs_microservice.models.nolabs_api_models_small_molecules_design_get_experiment_response import NolabsApiModelsSmallMoleculesDesignGetExperimentResponse

# TODO update the JSON string below
json = "{}"
# create an instance of NolabsApiModelsSmallMoleculesDesignGetExperimentResponse from a JSON string
nolabs_api_models_small_molecules_design_get_experiment_response_instance = NolabsApiModelsSmallMoleculesDesignGetExperimentResponse.from_json(json)
# print the JSON string representation of the object
print NolabsApiModelsSmallMoleculesDesignGetExperimentResponse.to_json()

# convert the object into a dict
nolabs_api_models_small_molecules_design_get_experiment_response_dict = nolabs_api_models_small_molecules_design_get_experiment_response_instance.to_dict()
# create an instance of NolabsApiModelsSmallMoleculesDesignGetExperimentResponse from a dict
nolabs_api_models_small_molecules_design_get_experiment_response_form_dict = nolabs_api_models_small_molecules_design_get_experiment_response.from_dict(nolabs_api_models_small_molecules_design_get_experiment_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


