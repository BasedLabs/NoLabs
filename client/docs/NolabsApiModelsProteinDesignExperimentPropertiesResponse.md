# NolabsApiModelsProteinDesignExperimentPropertiesResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pdb_file** | **object** |  | 
**pdb_file_name** | **object** |  | 
**contig** | **object** |  | 
**number_of_designs** | **object** |  | 
**hotspots** | [**Hotspots**](Hotspots.md) |  | [optional] 
**timesteps** | [**Timesteps**](Timesteps.md) |  | [optional] 

## Example

```python
from nolabs_microservice.models.nolabs_api_models_protein_design_experiment_properties_response import NolabsApiModelsProteinDesignExperimentPropertiesResponse

# TODO update the JSON string below
json = "{}"
# create an instance of NolabsApiModelsProteinDesignExperimentPropertiesResponse from a JSON string
nolabs_api_models_protein_design_experiment_properties_response_instance = NolabsApiModelsProteinDesignExperimentPropertiesResponse.from_json(json)
# print the JSON string representation of the object
print NolabsApiModelsProteinDesignExperimentPropertiesResponse.to_json()

# convert the object into a dict
nolabs_api_models_protein_design_experiment_properties_response_dict = nolabs_api_models_protein_design_experiment_properties_response_instance.to_dict()
# create an instance of NolabsApiModelsProteinDesignExperimentPropertiesResponse from a dict
nolabs_api_models_protein_design_experiment_properties_response_form_dict = nolabs_api_models_protein_design_experiment_properties_response.from_dict(nolabs_api_models_protein_design_experiment_properties_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


