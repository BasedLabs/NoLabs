# Connection


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**source_path** | **List[str]** |  | 
**target_path** | **List[str]** |  | 
**source_component_id** | **str** |  | 

## Example

```python
from biobuddy_microservice.models.connection import Connection

# TODO update the JSON string below
json = "{}"
# create an instance of Connection from a JSON string
connection_instance = Connection.from_json(json)
# print the JSON string representation of the object
print Connection.to_json()

# convert the object into a dict
connection_dict = connection_instance.to_dict()
# create an instance of Connection from a dict
connection_form_dict = connection.from_dict(connection_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


