# RunUmolPredictionRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**protein_sequence** | **str** |  | 
**ligand_smiles** | **str** |  | 
**msa_content** | **str** |  | 
**pocket_ids** | **List[int]** |  | 
**job_id** | **str** |  | [optional] 

## Example

```python
from umol_microservice.models.run_umol_prediction_request import RunUmolPredictionRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunUmolPredictionRequest from a JSON string
run_umol_prediction_request_instance = RunUmolPredictionRequest.from_json(json)
# print the JSON string representation of the object
print RunUmolPredictionRequest.to_json()

# convert the object into a dict
run_umol_prediction_request_dict = run_umol_prediction_request_instance.to_dict()
# create an instance of RunUmolPredictionRequest from a dict
run_umol_prediction_request_form_dict = run_umol_prediction_request.from_dict(run_umol_prediction_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


