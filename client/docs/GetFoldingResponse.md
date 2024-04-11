# GetFoldingResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pdb_contents** | [**PdbContents**](PdbContents.md) |  | 

## Example

```python
from nolabs_microservice.models.get_folding_response import GetFoldingResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetFoldingResponse from a JSON string
get_folding_response_instance = GetFoldingResponse.from_json(json)
# print the JSON string representation of the object
print GetFoldingResponse.to_json()

# convert the object into a dict
get_folding_response_dict = get_folding_response_instance.to_dict()
# create an instance of GetFoldingResponse from a dict
get_folding_response_form_dict = get_folding_response.from_dict(get_folding_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


