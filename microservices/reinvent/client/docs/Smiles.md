# Smiles


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**smiles** | **str** |  | 
**drug_likeness** | **float** |  | 
**score** | **float** |  | 
**stage** | **str** |  | 

## Example

```python
from reinvent_microservice.models.smiles import Smiles

# TODO update the JSON string below
json = "{}"
# create an instance of Smiles from a JSON string
smiles_instance = Smiles.from_json(json)
# print the JSON string representation of the object
print(Smiles.to_json())

# convert the object into a dict
smiles_dict = smiles_instance.to_dict()
# create an instance of Smiles from a dict
smiles_form_dict = smiles.from_dict(smiles_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


