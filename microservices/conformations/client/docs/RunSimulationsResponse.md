# RunSimulationsResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**errors** | **List[str]** |  | 
**pdb_content** | **str** |  | [optional] 

## Example

```python
from conformations_microservice.models.run_simulations_response import RunSimulationsResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunSimulationsResponse from a JSON string
run_simulations_response_instance = RunSimulationsResponse.from_json(json)
# print the JSON string representation of the object
print RunSimulationsResponse.to_json()

# convert the object into a dict
run_simulations_response_dict = run_simulations_response_instance.to_dict()
# create an instance of RunSimulationsResponse from a dict
run_simulations_response_form_dict = run_simulations_response.from_dict(run_simulations_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


