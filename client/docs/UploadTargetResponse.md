# UploadTargetResponse

Returns the list of targets since an uploaded .fasta file could contain multiple sequences inside

## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**result** | **object** |  | 

## Example

```python
from nolabs_microservice.models.upload_target_response import UploadTargetResponse

# TODO update the JSON string below
json = "{}"
# create an instance of UploadTargetResponse from a JSON string
upload_target_response_instance = UploadTargetResponse.from_json(json)
# print the JSON string representation of the object
print UploadTargetResponse.to_json()

# convert the object into a dict
upload_target_response_dict = upload_target_response_instance.to_dict()
# create an instance of UploadTargetResponse from a dict
upload_target_response_form_dict = upload_target_response.from_dict(upload_target_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


