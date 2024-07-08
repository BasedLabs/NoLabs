# WorkflowComponent


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**id** | **str** |  | 
**name** | **str** |  | 
**description** | **str** |  | 
**connections** | [**List[Connection]**](Connection.md) |  | 

## Example

```python
from biobuddy_microservice.models.workflow_component import WorkflowComponent

# TODO update the JSON string below
json = "{}"
# create an instance of WorkflowComponent from a JSON string
workflow_component_instance = WorkflowComponent.from_json(json)
# print the JSON string representation of the object
print WorkflowComponent.to_json()

# convert the object into a dict
workflow_component_dict = workflow_component_instance.to_dict()
# create an instance of WorkflowComponent from a dict
workflow_component_form_dict = workflow_component.from_dict(workflow_component_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


