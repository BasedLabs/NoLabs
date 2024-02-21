# SDFResult


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**sdf_file_name** | **str** |  | 
**sdf_content** | **str** |  | 
**confidence** | **float** |  | 
**scored_affinity** | **float** |  | 
**minimized_affinity** | **float** |  | 

## Example

```python
from diffdock_microservice.models.sdf_result import SDFResult

# TODO update the JSON string below
json = "{}"
# create an instance of SDFResult from a JSON string
sdf_result_instance = SDFResult.from_json(json)
# print the JSON string representation of the object
print SDFResult.to_json()

# convert the object into a dict
sdf_result_dict = sdf_result_instance.to_dict()
# create an instance of SDFResult from a dict
sdf_result_form_dict = sdf_result.from_dict(sdf_result_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


