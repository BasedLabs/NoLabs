# SamplingSizeRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**number_of_molecules_to_design** | **int** |  | 

## Example

```python
from reinvent_microservice.models.sampling_size_request import SamplingSizeRequest

# TODO update the JSON string below
json = "{}"
# create an instance of SamplingSizeRequest from a JSON string
sampling_size_request_instance = SamplingSizeRequest.from_json(json)
# print the JSON string representation of the object
print(SamplingSizeRequest.to_json())

# convert the object into a dict
sampling_size_request_dict = sampling_size_request_instance.to_dict()
# create an instance of SamplingSizeRequest from a dict
sampling_size_request_form_dict = sampling_size_request.from_dict(sampling_size_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


