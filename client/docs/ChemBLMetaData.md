# ChemBLMetaData


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**chembl_id** | **object** |  | 
**link** | **object** |  | 
**pref_name** | **object** |  | 

## Example

```python
from nolabs_microservice.models.chem_bl_meta_data import ChemBLMetaData

# TODO update the JSON string below
json = "{}"
# create an instance of ChemBLMetaData from a JSON string
chem_bl_meta_data_instance = ChemBLMetaData.from_json(json)
# print the JSON string representation of the object
print ChemBLMetaData.to_json()

# convert the object into a dict
chem_bl_meta_data_dict = chem_bl_meta_data_instance.to_dict()
# create an instance of ChemBLMetaData from a dict
chem_bl_meta_data_form_dict = chem_bl_meta_data.from_dict(chem_bl_meta_data_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


