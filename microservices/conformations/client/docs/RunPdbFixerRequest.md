# RunPdbFixerRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**job_id** | **str** |  | 
**pdb_content** | **str** |  | 
**replace_nonstandard_residues** | **bool** |  | [optional] [default to False]
**add_missing_atoms** | **bool** |  | [optional] [default to False]
**add_missing_hydrogens** | **bool** |  | [optional] [default to True]

## Example

```python
from conformations_microservice.models.run_pdb_fixer_request import RunPdbFixerRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunPdbFixerRequest from a JSON string
run_pdb_fixer_request_instance = RunPdbFixerRequest.from_json(json)
# print the JSON string representation of the object
print RunPdbFixerRequest.to_json()

# convert the object into a dict
run_pdb_fixer_request_dict = run_pdb_fixer_request_instance.to_dict()
# create an instance of RunPdbFixerRequest from a dict
run_pdb_fixer_request_form_dict = run_pdb_fixer_request.from_dict(run_pdb_fixer_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


