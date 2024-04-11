# ChemBLData


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**content** | [**Content**](Content.md) |  | [optional] 
**metadata** | [**ChemBLMetaData**](ChemBLMetaData.md) |  | 

## Example

```python
from nolabs_microservice.models.chem_bl_data import ChemBLData

# TODO update the JSON string below
json = "{}"
# create an instance of ChemBLData from a JSON string
chem_bl_data_instance = ChemBLData.from_json(json)
# print the JSON string representation of the object
print ChemBLData.to_json()

# convert the object into a dict
chem_bl_data_dict = chem_bl_data_instance.to_dict()
# create an instance of ChemBLData from a dict
chem_bl_data_form_dict = chem_bl_data.from_dict(chem_bl_data_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


