# Molecule


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**chembl_id** | **str** |  | 
**molecule_type** | **str** |  | 
**synonyms** | **List[str]** |  | 
**smiles** | **str** |  | 
**link** | **str** |  | 
**pref_name** | **str** |  | [optional] 

## Example

```python
from chembl_query_microservice.models.molecule import Molecule

# TODO update the JSON string below
json = "{}"
# create an instance of Molecule from a JSON string
molecule_instance = Molecule.from_json(json)
# print the JSON string representation of the object
print Molecule.to_json()

# convert the object into a dict
molecule_dict = molecule_instance.to_dict()
# create an instance of Molecule from a dict
molecule_form_dict = molecule.from_dict(molecule_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


