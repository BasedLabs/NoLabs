# SendMessageResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**biobuddy_response** | [**Message**](Message.md) |  | 

## Example

```python
from nolabs_microservice.models.send_message_response import SendMessageResponse

# TODO update the JSON string below
json = "{}"
# create an instance of SendMessageResponse from a JSON string
send_message_response_instance = SendMessageResponse.from_json(json)
# print the JSON string representation of the object
print SendMessageResponse.to_json()

# convert the object into a dict
send_message_response_dict = send_message_response_instance.to_dict()
# create an instance of SendMessageResponse from a dict
send_message_response_form_dict = send_message_response.from_dict(send_message_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


