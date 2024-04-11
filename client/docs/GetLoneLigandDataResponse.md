# GetLoneLigandDataResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**ligand_id** | **object** |  | 
**ligand_name** | **object** |  | 
**ligand_sdf** | **object** |  | 
**ligand_smiles** | **object** |  | 

## Example

```python
from nolabs_microservice.models.get_lone_ligand_data_response import GetLoneLigandDataResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetLoneLigandDataResponse from a JSON string
get_lone_ligand_data_response_instance = GetLoneLigandDataResponse.from_json(json)
# print the JSON string representation of the object
print GetLoneLigandDataResponse.to_json()

# convert the object into a dict
get_lone_ligand_data_response_dict = get_lone_ligand_data_response_instance.to_dict()
# create an instance of GetLoneLigandDataResponse from a dict
get_lone_ligand_data_response_form_dict = get_lone_ligand_data_response.from_dict(get_lone_ligand_data_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


