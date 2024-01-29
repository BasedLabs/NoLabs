# RunP2RankPredictionRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pdb_contents** | **str** |  | 
**job_id** | **str** |  | [optional] 

## Example

```python
from p2rank_microservice.models.run_p2_rank_prediction_request import RunP2RankPredictionRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunP2RankPredictionRequest from a JSON string
run_p2_rank_prediction_request_instance = RunP2RankPredictionRequest.from_json(json)
# print the JSON string representation of the object
print RunP2RankPredictionRequest.to_json()

# convert the object into a dict
run_p2_rank_prediction_request_dict = run_p2_rank_prediction_request_instance.to_dict()
# create an instance of RunP2RankPredictionRequest from a dict
run_p2_rank_prediction_request_form_dict = run_p2_rank_prediction_request.from_dict(run_p2_rank_prediction_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


