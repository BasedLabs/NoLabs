# GetTargetBindingPocketResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pocket_ids** | [**PocketIds**](PocketIds.md) |  | 

## Example

```python
from nolabs_microservice.models.get_target_binding_pocket_response import GetTargetBindingPocketResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetTargetBindingPocketResponse from a JSON string
get_target_binding_pocket_response_instance = GetTargetBindingPocketResponse.from_json(json)
# print the JSON string representation of the object
print GetTargetBindingPocketResponse.to_json()

# convert the object into a dict
get_target_binding_pocket_response_dict = get_target_binding_pocket_response_instance.to_dict()
# create an instance of GetTargetBindingPocketResponse from a dict
get_target_binding_pocket_response_form_dict = get_target_binding_pocket_response.from_dict(get_target_binding_pocket_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


