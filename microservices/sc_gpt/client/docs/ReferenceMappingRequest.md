# ReferenceMappingRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**query_path** | **str** |  | 
**reference_path** | **str** |  | 
**cell_type_key** | **str** |  | [optional] 
**gene_col** | **str** |  | [optional] [default to 'gene_name']
**batch_size** | **int** |  | [optional] [default to 64]
**device** | **str** |  | [optional] [default to 'cpu']
**k_neighbors** | **int** |  | [optional] [default to 10]
**calculate_metrics** | **bool** |  | [optional] [default to False]

## Example

```python
from sc_gpt_microservice.models.reference_mapping_request import ReferenceMappingRequest

# TODO update the JSON string below
json = "{}"
# create an instance of ReferenceMappingRequest from a JSON string
reference_mapping_request_instance = ReferenceMappingRequest.from_json(json)
# print the JSON string representation of the object
print(ReferenceMappingRequest.to_json())

# convert the object into a dict
reference_mapping_request_dict = reference_mapping_request_instance.to_dict()
# create an instance of ReferenceMappingRequest from a dict
reference_mapping_request_from_dict = ReferenceMappingRequest.from_dict(reference_mapping_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


