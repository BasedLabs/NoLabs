# SequenceQueryRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**sequence** | **str** |  | 
**sequence_type** | **str** |  | 
**identity_cutoff** | [**IdentityCutoff**](IdentityCutoff.md) |  | [optional] 
**evalue_cutoff** | [**EvalueCutoff**](EvalueCutoff.md) |  | [optional] 
**job_id** | [**JobId**](JobId.md) |  | [optional] 

## Example

```python
from rcsb_pdb_query_microservice.models.sequence_query_request import SequenceQueryRequest

# TODO update the JSON string below
json = "{}"
# create an instance of SequenceQueryRequest from a JSON string
sequence_query_request_instance = SequenceQueryRequest.from_json(json)
# print the JSON string representation of the object
print SequenceQueryRequest.to_json()

# convert the object into a dict
sequence_query_request_dict = sequence_query_request_instance.to_dict()
# create an instance of SequenceQueryRequest from a dict
sequence_query_request_form_dict = sequence_query_request.from_dict(sequence_query_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


