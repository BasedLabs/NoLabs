# PredictBindingPocketResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pocket_ids** | [**PocketIds**](PocketIds.md) |  | 

## Example

```python
from nolabs_microservice.models.predict_binding_pocket_response import PredictBindingPocketResponse

# TODO update the JSON string below
json = "{}"
# create an instance of PredictBindingPocketResponse from a JSON string
predict_binding_pocket_response_instance = PredictBindingPocketResponse.from_json(json)
# print the JSON string representation of the object
print PredictBindingPocketResponse.to_json()

# convert the object into a dict
predict_binding_pocket_response_dict = predict_binding_pocket_response_instance.to_dict()
# create an instance of PredictBindingPocketResponse from a dict
predict_binding_pocket_response_form_dict = predict_binding_pocket_response.from_dict(predict_binding_pocket_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


