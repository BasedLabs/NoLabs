# RunLocalisationPredictionResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**cytosolic_proteins** | **float** |  | 
**mitochondrial_proteins** | **float** |  | 
**nuclear_proteins** | **float** |  | 
**other_proteins** | **float** |  | 
**extracellular_secreted_proteins** | **float** |  | 
**errors** | **List[str]** |  | 

## Example

```python
from localisation_microservice.models.run_localisation_prediction_response import RunLocalisationPredictionResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunLocalisationPredictionResponse from a JSON string
run_localisation_prediction_response_instance = RunLocalisationPredictionResponse.from_json(json)
# print the JSON string representation of the object
print RunLocalisationPredictionResponse.to_json()

# convert the object into a dict
run_localisation_prediction_response_dict = run_localisation_prediction_response_instance.to_dict()
# create an instance of RunLocalisationPredictionResponse from a dict
run_localisation_prediction_response_form_dict = run_localisation_prediction_response.from_dict(run_localisation_prediction_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


