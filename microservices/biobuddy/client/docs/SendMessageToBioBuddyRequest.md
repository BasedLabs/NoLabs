# SendMessageToBioBuddyRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **str** |  | 
**message_content** | **str** |  | 
**previous_messages** | **List[Dict[str, str]]** |  | 
**available_components** | [**List[Component]**](Component.md) |  | 
**current_workflow** | [**List[WorkflowComponent]**](WorkflowComponent.md) |  | 
**tools** | **List[object]** |  | 
**job_id** | [**JobId**](JobId.md) |  | [optional] 

## Example

```python
from biobuddy_microservice.models.send_message_to_bio_buddy_request import SendMessageToBioBuddyRequest

# TODO update the JSON string below
json = "{}"
# create an instance of SendMessageToBioBuddyRequest from a JSON string
send_message_to_bio_buddy_request_instance = SendMessageToBioBuddyRequest.from_json(json)
# print the JSON string representation of the object
print SendMessageToBioBuddyRequest.to_json()

# convert the object into a dict
send_message_to_bio_buddy_request_dict = send_message_to_bio_buddy_request_instance.to_dict()
# create an instance of SendMessageToBioBuddyRequest from a dict
send_message_to_bio_buddy_request_form_dict = send_message_to_bio_buddy_request.from_dict(send_message_to_bio_buddy_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


