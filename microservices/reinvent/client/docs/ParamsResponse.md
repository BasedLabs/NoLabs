# ParamsResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**center_x** | **float** |  | 
**center_y** | **float** |  | 
**center_z** | **float** |  | 
**size_x** | **float** |  | 
**size_y** | **float** |  | 
**size_z** | **float** |  | 
**batch_size** | **float** |  | 
**minscore** | **float** |  | 
**epochs** | **float** |  | 

## Example

```python
from reinvent_microservice.models.params_response import ParamsResponse

# TODO update the JSON string below
json = "{}"
# create an instance of ParamsResponse from a JSON string
params_response_instance = ParamsResponse.from_json(json)
# print the JSON string representation of the object
print(ParamsResponse.to_json())

# convert the object into a dict
params_response_dict = params_response_instance.to_dict()
# create an instance of ParamsResponse from a dict
params_response_form_dict = params_response.from_dict(params_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


