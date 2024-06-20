# SequenceQuery

A query for a BLAST search.  - sequence: could be nucleotide sequence for blastn, tblastx, or tblastn, or amino acid sequence for blastp or blastx.

## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**sequence** | **str** |  | 
**type** | [**BlastType**](BlastType.md) |  | 
**descriptions** | **int** |  | [optional] 
**alignments** | **int** |  | [optional] 
**hitlist_size** | **int** |  | [optional] 
**expect** | **float** |  | [optional] 

## Example

```python
from blast_query_microservice.models.sequence_query import SequenceQuery

# TODO update the JSON string below
json = "{}"
# create an instance of SequenceQuery from a JSON string
sequence_query_instance = SequenceQuery.from_json(json)
# print the JSON string representation of the object
print(SequenceQuery.to_json())

# convert the object into a dict
sequence_query_dict = sequence_query_instance.to_dict()
# create an instance of SequenceQuery from a dict
sequence_query_from_dict = SequenceQuery.from_dict(sequence_query_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


