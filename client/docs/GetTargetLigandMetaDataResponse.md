# GetTargetLigandMetaDataResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**ligand_id** | **object** |  | 
**ligand_name** | **object** |  | 

## Example

```python
from nolabs_microservice.models.get_target_ligand_meta_data_response import GetTargetLigandMetaDataResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetTargetLigandMetaDataResponse from a JSON string
get_target_ligand_meta_data_response_instance = GetTargetLigandMetaDataResponse.from_json(json)
# print the JSON string representation of the object
print GetTargetLigandMetaDataResponse.to_json()

# convert the object into a dict
get_target_ligand_meta_data_response_dict = get_target_ligand_meta_data_response_instance.to_dict()
# create an instance of GetTargetLigandMetaDataResponse from a dict
get_target_ligand_meta_data_response_form_dict = get_target_ligand_meta_data_response.from_dict(get_target_ligand_meta_data_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


