# RunMsaPredictionRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**api_url** | **str** |  | 
**fasta_contents** | **str** |  | 
**job_id** | **str** |  | [optional] 

## Example

```python
from msa_light_microservice.models.run_msa_prediction_request import RunMsaPredictionRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunMsaPredictionRequest from a JSON string
run_msa_prediction_request_instance = RunMsaPredictionRequest.from_json(json)
# print the JSON string representation of the object
print RunMsaPredictionRequest.to_json()

# convert the object into a dict
run_msa_prediction_request_dict = run_msa_prediction_request_instance.to_dict()
# create an instance of RunMsaPredictionRequest from a dict
run_msa_prediction_request_form_dict = run_msa_prediction_request.from_dict(run_msa_prediction_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


