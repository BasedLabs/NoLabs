# GetJobBindingPocketDataResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pocket_ids** | [**PocketIds**](PocketIds.md) |  | 

## Example

```python
from nolabs_microservice.models.get_job_binding_pocket_data_response import GetJobBindingPocketDataResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetJobBindingPocketDataResponse from a JSON string
get_job_binding_pocket_data_response_instance = GetJobBindingPocketDataResponse.from_json(json)
# print the JSON string representation of the object
print GetJobBindingPocketDataResponse.to_json()

# convert the object into a dict
get_job_binding_pocket_data_response_dict = get_job_binding_pocket_data_response_instance.to_dict()
# create an instance of GetJobBindingPocketDataResponse from a dict
get_job_binding_pocket_data_response_form_dict = get_job_binding_pocket_data_response.from_dict(get_job_binding_pocket_data_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


