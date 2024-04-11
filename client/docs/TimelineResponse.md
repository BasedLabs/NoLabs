# TimelineResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**message** | **object** |  | 
**error** | [**Error**](Error.md) |  | 
**created_at** | **object** |  | 

## Example

```python
from nolabs_microservice.models.timeline_response import TimelineResponse

# TODO update the JSON string below
json = "{}"
# create an instance of TimelineResponse from a JSON string
timeline_response_instance = TimelineResponse.from_json(json)
# print the JSON string representation of the object
print TimelineResponse.to_json()

# convert the object into a dict
timeline_response_dict = timeline_response_instance.to_dict()
# create an instance of TimelineResponse from a dict
timeline_response_form_dict = timeline_response.from_dict(timeline_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


