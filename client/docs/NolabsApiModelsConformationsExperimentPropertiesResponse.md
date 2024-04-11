# NolabsApiModelsConformationsExperimentPropertiesResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pdb_file** | **object** |  | 
**pdb_file_name** | **object** |  | 
**total_frames** | **object** |  | 
**temperature_k** | **object** |  | 
**take_frame_every** | **object** |  | 
**step_size** | **object** |  | 
**replace_non_standard_residues** | **object** |  | 
**add_missing_atoms** | **object** |  | 
**add_missing_hydrogens** | **object** |  | 
**friction_coeff** | **object** |  | 
**ignore_missing_atoms** | **object** |  | 
**integrator** | [**IntegratorsRequest**](IntegratorsRequest.md) |  | [optional] 

## Example

```python
from nolabs_microservice.models.nolabs_api_models_conformations_experiment_properties_response import NolabsApiModelsConformationsExperimentPropertiesResponse

# TODO update the JSON string below
json = "{}"
# create an instance of NolabsApiModelsConformationsExperimentPropertiesResponse from a JSON string
nolabs_api_models_conformations_experiment_properties_response_instance = NolabsApiModelsConformationsExperimentPropertiesResponse.from_json(json)
# print the JSON string representation of the object
print NolabsApiModelsConformationsExperimentPropertiesResponse.to_json()

# convert the object into a dict
nolabs_api_models_conformations_experiment_properties_response_dict = nolabs_api_models_conformations_experiment_properties_response_instance.to_dict()
# create an instance of NolabsApiModelsConformationsExperimentPropertiesResponse from a dict
nolabs_api_models_conformations_experiment_properties_response_form_dict = nolabs_api_models_conformations_experiment_properties_response.from_dict(nolabs_api_models_conformations_experiment_properties_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


