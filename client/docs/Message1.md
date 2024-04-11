# Message1


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**content** | **object** |  | 
**function_name** | **object** |  | 
**parameters** | **object** |  | 
**data** | [**FunctionCallReturnData**](FunctionCallReturnData.md) |  | [optional] 

## Example

```python
from nolabs_microservice.models.message1 import Message1

# TODO update the JSON string below
json = "{}"
# create an instance of Message1 from a JSON string
message1_instance = Message1.from_json(json)
# print the JSON string representation of the object
print Message1.to_json()

# convert the object into a dict
message1_dict = message1_instance.to_dict()
# create an instance of Message1 from a dict
message1_form_dict = message1.from_dict(message1_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


