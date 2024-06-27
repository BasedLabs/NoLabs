# EmbedResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**embeddings** | **List[List[float]]** |  | 

## Example

```python
from sc_gpt_microservice.models.embed_response import EmbedResponse

# TODO update the JSON string below
json = "{}"
# create an instance of EmbedResponse from a JSON string
embed_response_instance = EmbedResponse.from_json(json)
# print the JSON string representation of the object
print(EmbedResponse.to_json())

# convert the object into a dict
embed_response_dict = embed_response_instance.to_dict()
# create an instance of EmbedResponse from a dict
embed_response_from_dict = EmbedResponse.from_dict(embed_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


