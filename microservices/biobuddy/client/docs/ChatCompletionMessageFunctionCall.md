# ChatCompletionMessageFunctionCall


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**arguments** | **object** |  | 
**name** | **object** |  | 

## Example

```python
from biobuddy_microservice.models.chat_completion_message_function_call import ChatCompletionMessageFunctionCall

# TODO update the JSON string below
json = "{}"
# create an instance of ChatCompletionMessageFunctionCall from a JSON string
chat_completion_message_function_call_instance = ChatCompletionMessageFunctionCall.from_json(json)
# print the JSON string representation of the object
print ChatCompletionMessageFunctionCall.to_json()

# convert the object into a dict
chat_completion_message_function_call_dict = chat_completion_message_function_call_instance.to_dict()
# create an instance of ChatCompletionMessageFunctionCall from a dict
chat_completion_message_function_call_form_dict = chat_completion_message_function_call.from_dict(chat_completion_message_function_call_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


