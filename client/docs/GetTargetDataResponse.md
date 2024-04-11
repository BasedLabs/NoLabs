# GetTargetDataResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**protein_sequence** | **object** |  | 
**protein_pdb** | [**ProteinPdb**](ProteinPdb.md) |  | [optional] 

## Example

```python
from nolabs_microservice.models.get_target_data_response import GetTargetDataResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetTargetDataResponse from a JSON string
get_target_data_response_instance = GetTargetDataResponse.from_json(json)
# print the JSON string representation of the object
print GetTargetDataResponse.to_json()

# convert the object into a dict
get_target_data_response_dict = get_target_data_response_instance.to_dict()
# create an instance of GetTargetDataResponse from a dict
get_target_data_response_form_dict = get_target_data_response.from_dict(get_target_data_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


