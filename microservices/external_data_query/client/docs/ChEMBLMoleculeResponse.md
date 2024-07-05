# ChEMBLMoleculeResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**molecules** | [**List[Molecule]**](Molecule.md) |  | 

## Example

```python
from external_data_query_microservice.models.ch_embl_molecule_response import ChEMBLMoleculeResponse

# TODO update the JSON string below
json = "{}"
# create an instance of ChEMBLMoleculeResponse from a JSON string
ch_embl_molecule_response_instance = ChEMBLMoleculeResponse.from_json(json)
# print the JSON string representation of the object
print ChEMBLMoleculeResponse.to_json()

# convert the object into a dict
ch_embl_molecule_response_dict = ch_embl_molecule_response_instance.to_dict()
# create an instance of ChEMBLMoleculeResponse from a dict
ch_embl_molecule_response_form_dict = ch_embl_molecule_response.from_dict(ch_embl_molecule_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


