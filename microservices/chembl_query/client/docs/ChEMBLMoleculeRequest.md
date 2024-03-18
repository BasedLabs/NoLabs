# ChEMBLMoleculeRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**filters** | **object** |  | [optional] 
**order_by** | **str** |  | [optional] 
**limit** | **int** |  | [optional] [default to 20]
**job_id** | [**JobId**](JobId.md) |  | [optional] 

## Example

```python
from chembl_query_microservice.models.ch_embl_molecule_request import ChEMBLMoleculeRequest

# TODO update the JSON string below
json = "{}"
# create an instance of ChEMBLMoleculeRequest from a JSON string
ch_embl_molecule_request_instance = ChEMBLMoleculeRequest.from_json(json)
# print the JSON string representation of the object
print ChEMBLMoleculeRequest.to_json()

# convert the object into a dict
ch_embl_molecule_request_dict = ch_embl_molecule_request_instance.to_dict()
# create an instance of ChEMBLMoleculeRequest from a dict
ch_embl_molecule_request_form_dict = ch_embl_molecule_request.from_dict(ch_embl_molecule_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


