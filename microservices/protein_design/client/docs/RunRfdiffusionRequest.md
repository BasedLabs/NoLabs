# RunRfdiffusionRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pdb_content** | **str** |  | 
**hotspots** | **str** |  | [optional] 
**contig** | **str** |  | [optional] 
**timesteps** | **int** |  | [optional] 
**number_of_designs** | **int** |  | [optional] 

## Example

```python
from protein_design_microservice.models.run_rfdiffusion_request import RunRfdiffusionRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunRfdiffusionRequest from a JSON string
run_rfdiffusion_request_instance = RunRfdiffusionRequest.from_json(json)
# print the JSON string representation of the object
print RunRfdiffusionRequest.to_json()

# convert the object into a dict
run_rfdiffusion_request_dict = run_rfdiffusion_request_instance.to_dict()
# create an instance of RunRfdiffusionRequest from a dict
run_rfdiffusion_request_form_dict = run_rfdiffusion_request.from_dict(run_rfdiffusion_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


