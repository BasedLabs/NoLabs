# GetFastaFilesBySearchQueryRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**search_query** | **str** |  | 
**max_results** | **int** |  | [optional] 
**exact_match** | **bool** |  | [optional] [default to False]
**job_id** | [**JobId**](JobId.md) |  | [optional] 

## Example

```python
from rcsb_pdb_query_microservice.models.get_fasta_files_by_search_query_request import GetFastaFilesBySearchQueryRequest

# TODO update the JSON string below
json = "{}"
# create an instance of GetFastaFilesBySearchQueryRequest from a JSON string
get_fasta_files_by_search_query_request_instance = GetFastaFilesBySearchQueryRequest.from_json(json)
# print the JSON string representation of the object
print GetFastaFilesBySearchQueryRequest.to_json()

# convert the object into a dict
get_fasta_files_by_search_query_request_dict = get_fasta_files_by_search_query_request_instance.to_dict()
# create an instance of GetFastaFilesBySearchQueryRequest from a dict
get_fasta_files_by_search_query_request_form_dict = get_fasta_files_by_search_query_request.from_dict(get_fasta_files_by_search_query_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


