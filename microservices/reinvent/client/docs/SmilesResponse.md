# SmilesResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**smiles** | [**List[Smiles]**](Smiles.md) |  | 

## Example

```python
from reinvent_microservice.models.smiles_response import SmilesResponse

# TODO update the JSON string below
json = "{}"
# create an instance of SmilesResponse from a JSON string
smiles_response_instance = SmilesResponse.from_json(json)
# print the JSON string representation of the object
print(SmilesResponse.to_json())

# convert the object into a dict
smiles_response_dict = smiles_response_instance.to_dict()
# create an instance of SmilesResponse from a dict
smiles_response_form_dict = smiles_response.from_dict(smiles_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


