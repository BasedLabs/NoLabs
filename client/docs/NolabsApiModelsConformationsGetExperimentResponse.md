# NolabsApiModelsConformationsGetExperimentResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **object** |  | 
**experiment_name** | **object** |  | 
**properties** | [**NolabsApiModelsConformationsExperimentPropertiesResponse**](NolabsApiModelsConformationsExperimentPropertiesResponse.md) |  | 
**timeline** | **object** |  | 
**pdb_file** | [**PdbFile**](PdbFile.md) |  | 

## Example

```python
from nolabs_microservice.models.nolabs_api_models_conformations_get_experiment_response import NolabsApiModelsConformationsGetExperimentResponse

# TODO update the JSON string below
json = "{}"
# create an instance of NolabsApiModelsConformationsGetExperimentResponse from a JSON string
nolabs_api_models_conformations_get_experiment_response_instance = NolabsApiModelsConformationsGetExperimentResponse.from_json(json)
# print the JSON string representation of the object
print NolabsApiModelsConformationsGetExperimentResponse.to_json()

# convert the object into a dict
nolabs_api_models_conformations_get_experiment_response_dict = nolabs_api_models_conformations_get_experiment_response_instance.to_dict()
# create an instance of NolabsApiModelsConformationsGetExperimentResponse from a dict
nolabs_api_models_conformations_get_experiment_response_form_dict = nolabs_api_models_conformations_get_experiment_response.from_dict(nolabs_api_models_conformations_get_experiment_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


