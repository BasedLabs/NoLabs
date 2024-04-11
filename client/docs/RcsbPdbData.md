# RcsbPdbData


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**content** | [**Content**](Content.md) |  | [optional] 
**metadata** | [**RcsbPdbMetaData**](RcsbPdbMetaData.md) |  | 

## Example

```python
from nolabs_microservice.models.rcsb_pdb_data import RcsbPdbData

# TODO update the JSON string below
json = "{}"
# create an instance of RcsbPdbData from a JSON string
rcsb_pdb_data_instance = RcsbPdbData.from_json(json)
# print the JSON string representation of the object
print RcsbPdbData.to_json()

# convert the object into a dict
rcsb_pdb_data_dict = rcsb_pdb_data_instance.to_dict()
# create an instance of RcsbPdbData from a dict
rcsb_pdb_data_form_dict = rcsb_pdb_data.from_dict(rcsb_pdb_data_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


