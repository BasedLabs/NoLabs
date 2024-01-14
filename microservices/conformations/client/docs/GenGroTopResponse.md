# GenGroTopResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**gro** | **str** |  | 
**top** | **str** |  | 
**errors** | **List[str]** |  | [optional] 

## Example

```python
from conformations_microservice.models.gen_gro_top_response import GenGroTopResponse

# TODO update the JSON string below
json = "{}"
# create an instance of GenGroTopResponse from a JSON string
gen_gro_top_response_instance = GenGroTopResponse.from_json(json)
# print the JSON string representation of the object
print GenGroTopResponse.to_json()

# convert the object into a dict
gen_gro_top_response_dict = gen_gro_top_response_instance.to_dict()
# create an instance of GenGroTopResponse from a dict
gen_gro_top_response_form_dict = gen_gro_top_response.from_dict(gen_gro_top_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


