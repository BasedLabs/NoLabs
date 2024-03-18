# FetchedProtein


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**fasta_contents** | **str** |  | 
**link** | **str** |  | 

## Example

```python
from rcsb_pdb_query_microservice.models.fetched_protein import FetchedProtein

# TODO update the JSON string below
json = "{}"
# create an instance of FetchedProtein from a JSON string
fetched_protein_instance = FetchedProtein.from_json(json)
# print the JSON string representation of the object
print FetchedProtein.to_json()

# convert the object into a dict
fetched_protein_dict = fetched_protein_instance.to_dict()
# create an instance of FetchedProtein from a dict
fetched_protein_form_dict = fetched_protein.from_dict(fetched_protein_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


