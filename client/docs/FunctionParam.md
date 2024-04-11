# FunctionParam


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**name** | **object** |  | 
**value** | **object** |  | [optional] 

## Example

```python
from nolabs_microservice.models.function_param import FunctionParam

# TODO update the JSON string below
json = "{}"
# create an instance of FunctionParam from a JSON string
function_param_instance = FunctionParam.from_json(json)
# print the JSON string representation of the object
print FunctionParam.to_json()

# convert the object into a dict
function_param_dict = function_param_instance.to_dict()
# create an instance of FunctionParam from a dict
function_param_form_dict = function_param.from_dict(function_param_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


