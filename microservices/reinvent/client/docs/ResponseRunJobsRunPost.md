# ResponseRunJobsRunPost


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**id** | **object** |  | 
**name** | **object** |  | 
**created_at** | **object** |  | 
**running** | **object** |  | 
**learning_completed** | **object** |  | 

## Example

```python
from reinvent_microservice.models.response_run_jobs_run_post import ResponseRunJobsRunPost

# TODO update the JSON string below
json = "{}"
# create an instance of ResponseRunJobsRunPost from a JSON string
response_run_jobs_run_post_instance = ResponseRunJobsRunPost.from_json(json)
# print the JSON string representation of the object
print(ResponseRunJobsRunPost.to_json())

# convert the object into a dict
response_run_jobs_run_post_dict = response_run_jobs_run_post_instance.to_dict()
# create an instance of ResponseRunJobsRunPost from a dict
response_run_jobs_run_post_form_dict = response_run_jobs_run_post.from_dict(response_run_jobs_run_post_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


