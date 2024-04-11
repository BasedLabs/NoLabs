# GetLoneLigandMetaDataResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**ligand_id** | **object** |  | 
**ligand_name** | **object** |  | 
**image** | [**Image**](Image.md) |  | [optional] 

## Example

```python
from nolabs_microservice.models.get_lone_ligand_meta_data_response import GetLoneLigandMetaDataResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetLoneLigandMetaDataResponse from a JSON string
get_lone_ligand_meta_data_response_instance = GetLoneLigandMetaDataResponse.from_json(json)
# print the JSON string representation of the object
print GetLoneLigandMetaDataResponse.to_json()

# convert the object into a dict
get_lone_ligand_meta_data_response_dict = get_lone_ligand_meta_data_response_instance.to_dict()
# create an instance of GetLoneLigandMetaDataResponse from a dict
get_lone_ligand_meta_data_response_form_dict = get_lone_ligand_meta_data_response.from_dict(get_lone_ligand_meta_data_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


