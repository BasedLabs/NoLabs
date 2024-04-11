# RunFoldingResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **object** |  | 
**experiment_name** | **object** |  | 
**amino_acids** | **object** |  | 

## Example

```python
from nolabs_microservice.models.run_folding_response import RunFoldingResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunFoldingResponse from a JSON string
run_folding_response_instance = RunFoldingResponse.from_json(json)
# print the JSON string representation of the object
print RunFoldingResponse.to_json()

# convert the object into a dict
run_folding_response_dict = run_folding_response_instance.to_dict()
# create an instance of RunFoldingResponse from a dict
run_folding_response_form_dict = run_folding_response.from_dict(run_folding_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


