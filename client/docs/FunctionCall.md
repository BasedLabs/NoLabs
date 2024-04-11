# FunctionCall


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**function_name** | **object** |  | 
**parameters** | **object** |  | 
**data** | [**FunctionCallReturnData**](FunctionCallReturnData.md) |  | [optional] 

## Example

```python
from nolabs_microservice.models.function_call import FunctionCall

# TODO update the JSON string below
json = "{}"
# create an instance of FunctionCall from a JSON string
function_call_instance = FunctionCall.from_json(json)
# print the JSON string representation of the object
print FunctionCall.to_json()

# convert the object into a dict
function_call_dict = function_call_instance.to_dict()
# create an instance of FunctionCall from a dict
function_call_form_dict = function_call.from_dict(function_call_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


