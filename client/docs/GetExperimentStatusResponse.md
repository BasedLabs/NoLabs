# GetExperimentStatusResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**running** | **object** |  | 
**sampling_allowed** | **object** |  | 

## Example

```python
from nolabs_microservice.models.get_experiment_status_response import GetExperimentStatusResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetExperimentStatusResponse from a JSON string
get_experiment_status_response_instance = GetExperimentStatusResponse.from_json(json)
# print the JSON string representation of the object
print GetExperimentStatusResponse.to_json()

# convert the object into a dict
get_experiment_status_response_dict = get_experiment_status_response_instance.to_dict()
# create an instance of GetExperimentStatusResponse from a dict
get_experiment_status_response_form_dict = get_experiment_status_response.from_dict(get_experiment_status_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


