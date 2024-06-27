# ReferenceMappingResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**predictions** | **List[str]** |  | 
**accuracy** | **float** |  | 
**precision** | **float** |  | 
**recall** | **float** |  | 
**macro_f1** | **float** |  | 

## Example

```python
from sc_gpt_microservice.models.reference_mapping_response import ReferenceMappingResponse

# TODO update the JSON string below
json = "{}"
# create an instance of ReferenceMappingResponse from a JSON string
reference_mapping_response_instance = ReferenceMappingResponse.from_json(json)
# print the JSON string representation of the object
print(ReferenceMappingResponse.to_json())

# convert the object into a dict
reference_mapping_response_dict = reference_mapping_response_instance.to_dict()
# create an instance of ReferenceMappingResponse from a dict
reference_mapping_response_from_dict = ReferenceMappingResponse.from_dict(reference_mapping_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


