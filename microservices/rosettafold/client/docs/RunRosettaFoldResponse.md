# RunRosettaFoldResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**pdb_content** | [**PdbContent**](PdbContent.md) |  | 
**errors** | **List[str]** |  | 

## Example

```python
from rosettafold_microservice.models.run_rosetta_fold_response import RunRosettaFoldResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunRosettaFoldResponse from a JSON string
run_rosetta_fold_response_instance = RunRosettaFoldResponse.from_json(json)
# print the JSON string representation of the object
print RunRosettaFoldResponse.to_json()

# convert the object into a dict
run_rosetta_fold_response_dict = run_rosetta_fold_response_instance.to_dict()
# create an instance of RunRosettaFoldResponse from a dict
run_rosetta_fold_response_form_dict = run_rosetta_fold_response.from_dict(run_rosetta_fold_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


