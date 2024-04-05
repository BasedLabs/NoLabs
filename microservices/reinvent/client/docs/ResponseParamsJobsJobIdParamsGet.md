# ResponseParamsJobsJobIdParamsGet


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**center_x** | **object** |  | 
**center_y** | **object** |  | 
**center_z** | **object** |  | 
**size_x** | **object** |  | 
**size_y** | **object** |  | 
**size_z** | **object** |  | 
**batch_size** | **object** |  | 
**minscore** | **object** |  | 
**epochs** | **object** |  | 

## Example

```python
from reinvent_microservice.models.response_params_jobs_job_id_params_get import ResponseParamsJobsJobIdParamsGet

# TODO update the JSON string below
json = "{}"
# create an instance of ResponseParamsJobsJobIdParamsGet from a JSON string
response_params_jobs_job_id_params_get_instance = ResponseParamsJobsJobIdParamsGet.from_json(json)
# print the JSON string representation of the object
print(ResponseParamsJobsJobIdParamsGet.to_json())

# convert the object into a dict
response_params_jobs_job_id_params_get_dict = response_params_jobs_job_id_params_get_instance.to_dict()
# create an instance of ResponseParamsJobsJobIdParamsGet from a dict
response_params_jobs_job_id_params_get_form_dict = response_params_jobs_job_id_params_get.from_dict(response_params_jobs_job_id_params_get_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


