# RunSolubilityResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **object** |  | 
**experiment_name** | **object** |  | 
**amino_acids** | **object** |  | 
**errors** | **object** |  | [optional] 

## Example

```python
from nolabs_microservice.models.run_solubility_response import RunSolubilityResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunSolubilityResponse from a JSON string
run_solubility_response_instance = RunSolubilityResponse.from_json(json)
# print the JSON string representation of the object
print RunSolubilityResponse.to_json()

# convert the object into a dict
run_solubility_response_dict = run_solubility_response_instance.to_dict()
# create an instance of RunSolubilityResponse from a dict
run_solubility_response_form_dict = run_solubility_response.from_dict(run_solubility_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


