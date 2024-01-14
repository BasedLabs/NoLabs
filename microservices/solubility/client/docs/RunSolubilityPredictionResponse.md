# RunSolubilityPredictionResponse

RunSolubilityPredictionResponse(errors: List[str], soluble_probability: Optional[float])

## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**errors** | **List[str]** |  | 
**soluble_probability** | **float** |  | 

## Example

```python
from solubility_microservice.models.run_solubility_prediction_response import RunSolubilityPredictionResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunSolubilityPredictionResponse from a JSON string
run_solubility_prediction_response_instance = RunSolubilityPredictionResponse.from_json(json)
# print the JSON string representation of the object
print RunSolubilityPredictionResponse.to_json()

# convert the object into a dict
run_solubility_prediction_response_dict = run_solubility_prediction_response_instance.to_dict()
# create an instance of RunSolubilityPredictionResponse from a dict
run_solubility_prediction_response_form_dict = run_solubility_prediction_response.from_dict(run_solubility_prediction_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


