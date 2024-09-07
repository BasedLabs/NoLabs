# RunEsmFoldPredictionRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**protein_sequence** | **str** |  | 
**job_id** | **str** |  | [optional] 

## Example

```python
from esmfold_light_microservice.models.run_esm_fold_prediction_request import RunEsmFoldPredictionRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunEsmFoldPredictionRequest from a JSON string
run_esm_fold_prediction_request_instance = RunEsmFoldPredictionRequest.from_json(json)
# print the JSON string representation of the object
print RunEsmFoldPredictionRequest.to_json()

# convert the object into a dict
run_esm_fold_prediction_request_dict = run_esm_fold_prediction_request_instance.to_dict()
# create an instance of RunEsmFoldPredictionRequest from a dict
run_esm_fold_prediction_request_form_dict = run_esm_fold_prediction_request.from_dict(run_esm_fold_prediction_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


