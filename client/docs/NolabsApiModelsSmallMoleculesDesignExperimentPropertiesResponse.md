# NolabsApiModelsSmallMoleculesDesignExperimentPropertiesResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**center_x** | **object** |  | 
**center_y** | **object** |  | 
**center_z** | **object** |  | 
**size_x** | **object** |  | 
**size_y** | **object** |  | 
**size_z** | **object** |  | 
**batch_size** | **object** |  | 
**minscore** | **object** |  | 
**epochs** | **object** |  | 
**pdb_file** | [**PdbFile**](PdbFile.md) |  | 
**pdb_file_name** | [**PdbFileName**](PdbFileName.md) |  | 

## Example

```python
from nolabs_microservice.models.nolabs_api_models_small_molecules_design_experiment_properties_response import NolabsApiModelsSmallMoleculesDesignExperimentPropertiesResponse

# TODO update the JSON string below
json = "{}"
# create an instance of NolabsApiModelsSmallMoleculesDesignExperimentPropertiesResponse from a JSON string
nolabs_api_models_small_molecules_design_experiment_properties_response_instance = NolabsApiModelsSmallMoleculesDesignExperimentPropertiesResponse.from_json(json)
# print the JSON string representation of the object
print NolabsApiModelsSmallMoleculesDesignExperimentPropertiesResponse.to_json()

# convert the object into a dict
nolabs_api_models_small_molecules_design_experiment_properties_response_dict = nolabs_api_models_small_molecules_design_experiment_properties_response_instance.to_dict()
# create an instance of NolabsApiModelsSmallMoleculesDesignExperimentPropertiesResponse from a dict
nolabs_api_models_small_molecules_design_experiment_properties_response_form_dict = nolabs_api_models_small_molecules_design_experiment_properties_response.from_dict(nolabs_api_models_small_molecules_design_experiment_properties_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


