# RunRfdiffusionResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pdbs_content** | **List[str]** |  | [optional] 
**errors** | **List[str]** |  | [optional] 

## Example

```python
from protein_design_microservice.models.run_rfdiffusion_response import RunRfdiffusionResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunRfdiffusionResponse from a JSON string
run_rfdiffusion_response_instance = RunRfdiffusionResponse.from_json(json)
# print the JSON string representation of the object
print RunRfdiffusionResponse.to_json()

# convert the object into a dict
run_rfdiffusion_response_dict = run_rfdiffusion_response_instance.to_dict()
# create an instance of RunRfdiffusionResponse from a dict
run_rfdiffusion_response_form_dict = run_rfdiffusion_response.from_dict(run_rfdiffusion_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


