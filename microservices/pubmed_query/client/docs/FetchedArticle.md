# FetchedArticle


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**title** | **str** |  | 
**summary** | **str** |  | 
**link** | **str** |  | 

## Example

```python
from pubmed_query_microservice.models.fetched_article import FetchedArticle

# TODO update the JSON string below
json = "{}"
# create an instance of FetchedArticle from a JSON string
fetched_article_instance = FetchedArticle.from_json(json)
# print the JSON string representation of the object
print FetchedArticle.to_json()

# convert the object into a dict
fetched_article_dict = fetched_article_instance.to_dict()
# create an instance of FetchedArticle from a dict
fetched_article_form_dict = fetched_article.from_dict(fetched_article_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


