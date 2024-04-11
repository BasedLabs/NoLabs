# RunUmolDockingJobResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**predicted_pdb** | **object** |  | 
**predicted_sdf** | **object** |  | 
**plddt_array** | **object** |  | 
**job_id** | **object** |  | 

## Example

```python
from nolabs_microservice.models.run_umol_docking_job_response import RunUmolDockingJobResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunUmolDockingJobResponse from a JSON string
run_umol_docking_job_response_instance = RunUmolDockingJobResponse.from_json(json)
# print the JSON string representation of the object
print RunUmolDockingJobResponse.to_json()

# convert the object into a dict
run_umol_docking_job_response_dict = run_umol_docking_job_response_instance.to_dict()
# create an instance of RunUmolDockingJobResponse from a dict
run_umol_docking_job_response_form_dict = run_umol_docking_job_response.from_dict(run_umol_docking_job_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


