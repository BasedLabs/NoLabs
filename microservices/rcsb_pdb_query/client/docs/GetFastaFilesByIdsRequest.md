# GetFastaFilesByIdsRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**rcsb_pdb_ids** | **List[str]** |  | 
**job_id** | [**JobId**](JobId.md) |  | [optional] 

## Example

```python
from rcsb_pdb_query_microservice.models.get_fasta_files_by_ids_request import GetFastaFilesByIdsRequest

# TODO update the JSON string below
json = "{}"
# create an instance of GetFastaFilesByIdsRequest from a JSON string
get_fasta_files_by_ids_request_instance = GetFastaFilesByIdsRequest.from_json(json)
# print the JSON string representation of the object
print GetFastaFilesByIdsRequest.to_json()

# convert the object into a dict
get_fasta_files_by_ids_request_dict = get_fasta_files_by_ids_request_instance.to_dict()
# create an instance of GetFastaFilesByIdsRequest from a dict
get_fasta_files_by_ids_request_form_dict = get_fasta_files_by_ids_request.from_dict(get_fasta_files_by_ids_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


