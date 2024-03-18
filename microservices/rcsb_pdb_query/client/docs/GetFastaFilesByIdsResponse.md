# GetFastaFilesByIdsResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**fasta_contents** | [**List[FetchedProtein]**](FetchedProtein.md) |  | 

## Example

```python
from rcsb_pdb_query_microservice.models.get_fasta_files_by_ids_response import GetFastaFilesByIdsResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GetFastaFilesByIdsResponse from a JSON string
get_fasta_files_by_ids_response_instance = GetFastaFilesByIdsResponse.from_json(json)
# print the JSON string representation of the object
print GetFastaFilesByIdsResponse.to_json()

# convert the object into a dict
get_fasta_files_by_ids_response_dict = get_fasta_files_by_ids_response_instance.to_dict()
# create an instance of GetFastaFilesByIdsResponse from a dict
get_fasta_files_by_ids_response_form_dict = get_fasta_files_by_ids_response.from_dict(get_fasta_files_by_ids_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


