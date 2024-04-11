# PredictFoldingResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pdb_content** | [**PdbContent**](PdbContent.md) |  | [optional] 

## Example

```python
from nolabs_microservice.models.predict_folding_response import PredictFoldingResponse

# TODO update the JSON string below
json = "{}"
# create an instance of PredictFoldingResponse from a JSON string
predict_folding_response_instance = PredictFoldingResponse.from_json(json)
# print the JSON string representation of the object
print PredictFoldingResponse.to_json()

# convert the object into a dict
predict_folding_response_dict = predict_folding_response_instance.to_dict()
# create an instance of PredictFoldingResponse from a dict
predict_folding_response_form_dict = predict_folding_response.from_dict(predict_folding_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


