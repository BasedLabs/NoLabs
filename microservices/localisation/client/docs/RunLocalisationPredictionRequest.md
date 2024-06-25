# RunLocalisationPredictionRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**job_id** | **str** |  | 
**amino_acid_sequence** | **str** |  | 

## Example

```python
from localisation_microservice.models.run_localisation_prediction_request import RunLocalisationPredictionRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunLocalisationPredictionRequest from a JSON string
run_localisation_prediction_request_instance = RunLocalisationPredictionRequest.from_json(json)
# print the JSON string representation of the object
print RunLocalisationPredictionRequest.to_json()

# convert the object into a dict
run_localisation_prediction_request_dict = run_localisation_prediction_request_instance.to_dict()
# create an instance of RunLocalisationPredictionRequest from a dict
run_localisation_prediction_request_form_dict = run_localisation_prediction_request.from_dict(run_localisation_prediction_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


