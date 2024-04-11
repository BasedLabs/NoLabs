# RunDiffDockDockingJobResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**predicted_pdb** | **object** |  | 
**predicted_ligands** | **object** |  | 

## Example

```python
from nolabs_microservice.models.run_diff_dock_docking_job_response import RunDiffDockDockingJobResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunDiffDockDockingJobResponse from a JSON string
run_diff_dock_docking_job_response_instance = RunDiffDockDockingJobResponse.from_json(json)
# print the JSON string representation of the object
print RunDiffDockDockingJobResponse.to_json()

# convert the object into a dict
run_diff_dock_docking_job_response_dict = run_diff_dock_docking_job_response_instance.to_dict()
# create an instance of RunDiffDockDockingJobResponse from a dict
run_diff_dock_docking_job_response_form_dict = run_diff_dock_docking_job_response.from_dict(run_diff_dock_docking_job_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


