# ExperimentMetadataResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**id** | **object** |  | 
**name** | **object** |  | 
**var_date** | **object** |  | 

## Example

```python
from nolabs_microservice.models.experiment_metadata_response import ExperimentMetadataResponse

# TODO update the JSON string below
json = "{}"
# create an instance of ExperimentMetadataResponse from a JSON string
experiment_metadata_response_instance = ExperimentMetadataResponse.from_json(json)
# print the JSON string representation of the object
print ExperimentMetadataResponse.to_json()

# convert the object into a dict
experiment_metadata_response_dict = experiment_metadata_response_instance.to_dict()
# create an instance of ExperimentMetadataResponse from a dict
experiment_metadata_response_form_dict = experiment_metadata_response.from_dict(experiment_metadata_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


