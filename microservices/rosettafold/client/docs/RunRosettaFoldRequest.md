# RunRosettaFoldRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**fasta_content** | [**FastaContent**](FastaContent.md) |  | 
**a3m_content** | [**A3MContent**](A3MContent.md) |  | 

## Example

```python
from rosettafold_microservice.models.run_rosetta_fold_request import RunRosettaFoldRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunRosettaFoldRequest from a JSON string
run_rosetta_fold_request_instance = RunRosettaFoldRequest.from_json(json)
# print the JSON string representation of the object
print RunRosettaFoldRequest.to_json()

# convert the object into a dict
run_rosetta_fold_request_dict = run_rosetta_fold_request_instance.to_dict()
# create an instance of RunRosettaFoldRequest from a dict
run_rosetta_fold_request_form_dict = run_rosetta_fold_request.from_dict(run_rosetta_fold_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


