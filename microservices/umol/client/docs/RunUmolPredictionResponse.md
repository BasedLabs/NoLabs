# RunUmolPredictionResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**errors** | **List[str]** |  | 
**sdf_contents** | **str** |  | 
**plddt_array** | **List[int]** |  | 

## Example

```python
from umol_microservice.models.run_umol_prediction_response import RunUmolPredictionResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunUmolPredictionResponse from a JSON string
run_umol_prediction_response_instance = RunUmolPredictionResponse.from_json(json)
# print the JSON string representation of the object
print RunUmolPredictionResponse.to_json()

# convert the object into a dict
run_umol_prediction_response_dict = run_umol_prediction_response_instance.to_dict()
# create an instance of RunUmolPredictionResponse from a dict
run_umol_prediction_response_form_dict = run_umol_prediction_response.from_dict(run_umol_prediction_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


