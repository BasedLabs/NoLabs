# GoConfidenceResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**name** | **str** |  | 
**confidence** | **float** |  | 

## Example

```python
from gene_ontology_microservice.models.go_confidence_response import GoConfidenceResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GoConfidenceResponse from a JSON string
go_confidence_response_instance = GoConfidenceResponse.from_json(json)
# print the JSON string representation of the object
print GoConfidenceResponse.to_json()

# convert the object into a dict
go_confidence_response_dict = go_confidence_response_instance.to_dict()
# create an instance of GoConfidenceResponse from a dict
go_confidence_response_form_dict = go_confidence_response.from_dict(go_confidence_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


