# GetDockingParamsResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**folding_method** | **object** |  | 
**docking_method** | **object** |  | 

## Example

```python
from nolabs_microservice.models.get_docking_params_response import GetDockingParamsResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetDockingParamsResponse from a JSON string
get_docking_params_response_instance = GetDockingParamsResponse.from_json(json)
# print the JSON string representation of the object
print GetDockingParamsResponse.to_json()

# convert the object into a dict
get_docking_params_response_dict = get_docking_params_response_instance.to_dict()
# create an instance of GetDockingParamsResponse from a dict
get_docking_params_response_form_dict = get_docking_params_response.from_dict(get_docking_params_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


