# GetFastaFilesResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**fasta_contents** | [**List[FetchedProtein]**](FetchedProtein.md) |  | 

## Example

```python
from rcsb_pdb_query_microservice.models.get_fasta_files_response import GetFastaFilesResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetFastaFilesResponse from a JSON string
get_fasta_files_response_instance = GetFastaFilesResponse.from_json(json)
# print the JSON string representation of the object
print GetFastaFilesResponse.to_json()

# convert the object into a dict
get_fasta_files_response_dict = get_fasta_files_response_instance.to_dict()
# create an instance of GetFastaFilesResponse from a dict
get_fasta_files_response_form_dict = get_fasta_files_response.from_dict(get_fasta_files_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


