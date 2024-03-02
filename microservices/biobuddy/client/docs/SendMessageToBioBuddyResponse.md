# SendMessageToBioBuddyResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**chatgpt_reply** | [**ChatCompletionMessage**](ChatCompletionMessage.md) |  | 

## Example

```python
from biobuddy_microservice.models.send_message_to_bio_buddy_response import SendMessageToBioBuddyResponse

# TODO update the JSON string below
json = "{}"
# create an instance of SendMessageToBioBuddyResponse from a JSON string
send_message_to_bio_buddy_response_instance = SendMessageToBioBuddyResponse.from_json(json)
# print the JSON string representation of the object
print SendMessageToBioBuddyResponse.to_json()

# convert the object into a dict
send_message_to_bio_buddy_response_dict = send_message_to_bio_buddy_response_instance.to_dict()
# create an instance of SendMessageToBioBuddyResponse from a dict
send_message_to_bio_buddy_response_form_dict = send_message_to_bio_buddy_response.from_dict(send_message_to_bio_buddy_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


