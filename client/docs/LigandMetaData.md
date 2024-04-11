# LigandMetaData


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**ligand_id** | **object** |  | 
**ligand_name** | **object** |  | 
**description** | [**Description**](Description.md) |  | [optional] 
**link** | [**Link**](Link.md) |  | [optional] 
**image** | [**Image**](Image.md) |  | [optional] 

## Example

```python
from nolabs_microservice.models.ligand_meta_data import LigandMetaData

# TODO update the JSON string below
json = "{}"
# create an instance of LigandMetaData from a JSON string
ligand_meta_data_instance = LigandMetaData.from_json(json)
# print the JSON string representation of the object
print LigandMetaData.to_json()

# convert the object into a dict
ligand_meta_data_dict = ligand_meta_data_instance.to_dict()
# create an instance of LigandMetaData from a dict
ligand_meta_data_form_dict = ligand_meta_data.from_dict(ligand_meta_data_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


