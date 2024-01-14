# RunSolubilityPredictionRequest

RunSolubilityPredictionRequest(amino_acid_sequence: str)

## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**amino_acid_sequence** | **str** |  | 

## Example

```python
from solubility_microservice.models.run_solubility_prediction_request import RunSolubilityPredictionRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunSolubilityPredictionRequest from a JSON string
run_solubility_prediction_request_instance = RunSolubilityPredictionRequest.from_json(json)
# print the JSON string representation of the object
print RunSolubilityPredictionRequest.to_json()

# convert the object into a dict
run_solubility_prediction_request_dict = run_solubility_prediction_request_instance.to_dict()
# create an instance of RunSolubilityPredictionRequest from a dict
run_solubility_prediction_request_form_dict = run_solubility_prediction_request.from_dict(run_solubility_prediction_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


