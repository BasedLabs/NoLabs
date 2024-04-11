# RunProteinDesignResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **object** |  | 
**experiment_name** | **object** |  | 
**pdb_files** | **object** |  | 

## Example

```python
from nolabs_microservice.models.run_protein_design_response import RunProteinDesignResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunProteinDesignResponse from a JSON string
run_protein_design_response_instance = RunProteinDesignResponse.from_json(json)
# print the JSON string representation of the object
print RunProteinDesignResponse.to_json()

# convert the object into a dict
run_protein_design_response_dict = run_protein_design_response_instance.to_dict()
# create an instance of RunProteinDesignResponse from a dict
run_protein_design_response_form_dict = run_protein_design_response.from_dict(run_protein_design_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


