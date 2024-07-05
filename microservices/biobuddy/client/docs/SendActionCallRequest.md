# SendActionCallRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**action_text** | **str** |  | 
**tools** | **List[object]** |  | 
**job_id** | [**JobId**](JobId.md) |  | [optional] 

## Example

```python
from biobuddy_microservice.models.send_action_call_request import SendActionCallRequest

# TODO update the JSON string below
json = "{}"
# create an instance of SendActionCallRequest from a JSON string
send_action_call_request_instance = SendActionCallRequest.from_json(json)
# print the JSON string representation of the object
print SendActionCallRequest.to_json()

# convert the object into a dict
send_action_call_request_dict = send_action_call_request_instance.to_dict()
# create an instance of SendActionCallRequest from a dict
send_action_call_request_form_dict = send_action_call_request.from_dict(send_action_call_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


