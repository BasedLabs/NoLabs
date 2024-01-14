# GenGroTopRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**force_field** | [**GromacsForceFields**](GromacsForceFields.md) |  | 
**water_force_field** | [**GromacsWaterForceFields**](GromacsWaterForceFields.md) |  | 
**pdb_content** | **str** |  | 
**ignore_missing_atoms** | **bool** |  | [optional] [default to False]

## Example

```python
from conformations_microservice.models.gen_gro_top_request import GenGroTopRequest

# TODO update the JSON string below
json = "{}"
# create an instance of GenGroTopRequest from a JSON string
gen_gro_top_request_instance = GenGroTopRequest.from_json(json)
# print the JSON string representation of the object
print GenGroTopRequest.to_json()

# convert the object into a dict
gen_gro_top_request_dict = gen_gro_top_request_instance.to_dict()
# create an instance of GenGroTopRequest from a dict
gen_gro_top_request_form_dict = gen_gro_top_request.from_dict(gen_gro_top_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


