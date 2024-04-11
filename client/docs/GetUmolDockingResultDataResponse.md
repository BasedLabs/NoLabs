# GetUmolDockingResultDataResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**predicted_pdb** | **object** |  | 
**predicted_sdf** | **object** |  | 
**plddt_array** | **object** |  | 
**job_id** | **object** |  | 

## Example

```python
from nolabs_microservice.models.get_umol_docking_result_data_response import GetUmolDockingResultDataResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetUmolDockingResultDataResponse from a JSON string
get_umol_docking_result_data_response_instance = GetUmolDockingResultDataResponse.from_json(json)
# print the JSON string representation of the object
print GetUmolDockingResultDataResponse.to_json()

# convert the object into a dict
get_umol_docking_result_data_response_dict = get_umol_docking_result_data_response_instance.to_dict()
# create an instance of GetUmolDockingResultDataResponse from a dict
get_umol_docking_result_data_response_form_dict = get_umol_docking_result_data_response.from_dict(get_umol_docking_result_data_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


