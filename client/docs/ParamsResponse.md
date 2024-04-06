# ParamsResponse


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

## Example

```python
from reinvent_microservice.models.params_response import ParamsResponse

# TODO update the JSON string below
json = "{}"
# create an instance of ParamsResponse from a JSON string
params_response_instance = ParamsResponse.from_json(json)
# print the JSON string representation of the object
print ParamsResponse.to_json()

# convert the object into a dict
params_response_dict = params_response_instance.to_dict()
# create an instance of ParamsResponse from a dict
params_response_form_dict = params_response.from_dict(params_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


