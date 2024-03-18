# PubMedSearchResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**articles** | [**List[FetchedArticle]**](FetchedArticle.md) |  | 

## Example

```python
from pubmed_query_microservice.models.pub_med_search_response import PubMedSearchResponse

# TODO update the JSON string below
json = "{}"
# create an instance of PubMedSearchResponse from a JSON string
pub_med_search_response_instance = PubMedSearchResponse.from_json(json)
# print the JSON string representation of the object
print PubMedSearchResponse.to_json()

# convert the object into a dict
pub_med_search_response_dict = pub_med_search_response_instance.to_dict()
# create an instance of PubMedSearchResponse from a dict
pub_med_search_response_form_dict = pub_med_search_response.from_dict(pub_med_search_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


