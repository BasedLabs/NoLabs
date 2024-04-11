# GetDiffDockDockingResultDataResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**predicted_pdb** | **object** |  | 
**predicted_ligands** | **object** |  | 

## Example

```python
from nolabs_microservice.models.get_diff_dock_docking_result_data_response import GetDiffDockDockingResultDataResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetDiffDockDockingResultDataResponse from a JSON string
get_diff_dock_docking_result_data_response_instance = GetDiffDockDockingResultDataResponse.from_json(json)
# print the JSON string representation of the object
print GetDiffDockDockingResultDataResponse.to_json()

# convert the object into a dict
get_diff_dock_docking_result_data_response_dict = get_diff_dock_docking_result_data_response_instance.to_dict()
# create an instance of GetDiffDockDockingResultDataResponse from a dict
get_diff_dock_docking_result_data_response_form_dict = get_diff_dock_docking_result_data_response.from_dict(get_diff_dock_docking_result_data_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


