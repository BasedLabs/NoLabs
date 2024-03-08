# RunDiffDockPredictionRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pdb_contents** | **str** |  | 
**sdf_contents** | **str** |  | 
**inference_steps** | **int** |  | [optional] [default to 20]
**samples_per_complex** | **int** |  | [optional] [default to 40]
**batch_size** | **int** |  | [optional] [default to 10]
**actual_steps** | **int** |  | [optional] [default to 18]
**no_final_step_noise** | **bool** |  | [optional] [default to True]
**job_id** | **str** |  | [optional] 

## Example

```python
from diffdock_microservice.models.run_diff_dock_prediction_request import RunDiffDockPredictionRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunDiffDockPredictionRequest from a JSON string
run_diff_dock_prediction_request_instance = RunDiffDockPredictionRequest.from_json(json)
# print the JSON string representation of the object
print RunDiffDockPredictionRequest.to_json()

# convert the object into a dict
run_diff_dock_prediction_request_dict = run_diff_dock_prediction_request_instance.to_dict()
# create an instance of RunDiffDockPredictionRequest from a dict
run_diff_dock_prediction_request_form_dict = run_diff_dock_prediction_request.from_dict(run_diff_dock_prediction_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


