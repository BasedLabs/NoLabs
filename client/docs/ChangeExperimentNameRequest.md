# ChangeExperimentNameRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**id** | **object** |  | 
**name** | **object** |  | 

## Example

```python
from nolabs_microservice.models.change_experiment_name_request import ChangeExperimentNameRequest

# TODO update the JSON string below
json = "{}"
# create an instance of ChangeExperimentNameRequest from a JSON string
change_experiment_name_request_instance = ChangeExperimentNameRequest.from_json(json)
# print the JSON string representation of the object
print ChangeExperimentNameRequest.to_json()

# convert the object into a dict
change_experiment_name_request_dict = change_experiment_name_request_instance.to_dict()
# create an instance of ChangeExperimentNameRequest from a dict
change_experiment_name_request_form_dict = change_experiment_name_request.from_dict(change_experiment_name_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


