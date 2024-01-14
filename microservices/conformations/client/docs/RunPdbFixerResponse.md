# RunPdbFixerResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**errors** | **List[str]** |  | [optional] 
**pdb_content** | **str** |  | [optional] 

## Example

```python
from conformations_microservice.models.run_pdb_fixer_response import RunPdbFixerResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunPdbFixerResponse from a JSON string
run_pdb_fixer_response_instance = RunPdbFixerResponse.from_json(json)
# print the JSON string representation of the object
print RunPdbFixerResponse.to_json()

# convert the object into a dict
run_pdb_fixer_response_dict = run_pdb_fixer_response_instance.to_dict()
# create an instance of RunPdbFixerResponse from a dict
run_pdb_fixer_response_form_dict = run_pdb_fixer_response.from_dict(run_pdb_fixer_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


