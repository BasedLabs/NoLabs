# DiffDockLigandMetaData


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**job_id** | **object** |  | 
**target_id** | **object** |  | 
**ligand_id** | **object** |  | 
**predicted_ligand_file_name** | **object** |  | 
**minimized_affinity** | **object** |  | 
**scored_affinity** | **object** |  | 
**confidence** | [**Confidence**](Confidence.md) |  | [optional] 

## Example

```python
from nolabs_microservice.models.diff_dock_ligand_meta_data import DiffDockLigandMetaData

# TODO update the JSON string below
json = "{}"
# create an instance of DiffDockLigandMetaData from a JSON string
diff_dock_ligand_meta_data_instance = DiffDockLigandMetaData.from_json(json)
# print the JSON string representation of the object
print DiffDockLigandMetaData.to_json()

# convert the object into a dict
diff_dock_ligand_meta_data_dict = diff_dock_ligand_meta_data_instance.to_dict()
# create an instance of DiffDockLigandMetaData from a dict
diff_dock_ligand_meta_data_form_dict = diff_dock_ligand_meta_data.from_dict(diff_dock_ligand_meta_data_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


