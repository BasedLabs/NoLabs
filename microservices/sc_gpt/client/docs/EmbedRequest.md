# EmbedRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**dataset_path** | **str** |  | 
**gene_col** | **str** |  | [optional] [default to 'gene_name']
**batch_size** | **int** |  | [optional] [default to 64]
**device** | **str** |  | [optional] [default to 'cpu']

## Example

```python
from sc_gpt_microservice.models.embed_request import EmbedRequest

# TODO update the JSON string below
json = "{}"
# create an instance of EmbedRequest from a JSON string
embed_request_instance = EmbedRequest.from_json(json)
# print the JSON string representation of the object
print(EmbedRequest.to_json())

# convert the object into a dict
embed_request_dict = embed_request_instance.to_dict()
# create an instance of EmbedRequest from a dict
embed_request_from_dict = EmbedRequest.from_dict(embed_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


