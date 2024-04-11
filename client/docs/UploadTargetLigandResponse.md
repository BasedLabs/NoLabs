# UploadTargetLigandResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**ligand_meta_data** | [**LigandMetaData**](LigandMetaData.md) |  | 

## Example

```python
from nolabs_microservice.models.upload_target_ligand_response import UploadTargetLigandResponse

# TODO update the JSON string below
json = "{}"
# create an instance of UploadTargetLigandResponse from a JSON string
upload_target_ligand_response_instance = UploadTargetLigandResponse.from_json(json)
# print the JSON string representation of the object
print UploadTargetLigandResponse.to_json()

# convert the object into a dict
upload_target_ligand_response_dict = upload_target_ligand_response_instance.to_dict()
# create an instance of UploadTargetLigandResponse from a dict
upload_target_ligand_response_form_dict = upload_target_ligand_response.from_dict(upload_target_ligand_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


