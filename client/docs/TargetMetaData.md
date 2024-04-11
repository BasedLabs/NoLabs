# TargetMetaData


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**target_id** | **object** |  | 
**target_name** | **object** |  | 
**link** | [**Link**](Link.md) |  | [optional] 
**folding_method** | [**FoldingMethod**](FoldingMethod.md) |  | [optional] 

## Example

```python
from nolabs_microservice.models.target_meta_data import TargetMetaData

# TODO update the JSON string below
json = "{}"
# create an instance of TargetMetaData from a JSON string
target_meta_data_instance = TargetMetaData.from_json(json)
# print the JSON string representation of the object
print TargetMetaData.to_json()

# convert the object into a dict
target_meta_data_dict = target_meta_data_instance.to_dict()
# create an instance of TargetMetaData from a dict
target_meta_data_form_dict = target_meta_data.from_dict(target_meta_data_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


