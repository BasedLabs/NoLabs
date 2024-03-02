# RunDiffDockPredictionResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**success** | **bool** |  | 
**message** | **str** |  | 
**pdb_contents** | **str** |  | 
**sdf_results** | [**List[SDFResult]**](SDFResult.md) |  | 

## Example

```python
from diffdock_microservice.models.run_diff_dock_prediction_response import RunDiffDockPredictionResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunDiffDockPredictionResponse from a JSON string
run_diff_dock_prediction_response_instance = RunDiffDockPredictionResponse.from_json(json)
# print the JSON string representation of the object
print RunDiffDockPredictionResponse.to_json()

# convert the object into a dict
run_diff_dock_prediction_response_dict = run_diff_dock_prediction_response_instance.to_dict()
# create an instance of RunDiffDockPredictionResponse from a dict
run_diff_dock_prediction_response_form_dict = run_diff_dock_prediction_response.from_dict(run_diff_dock_prediction_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


