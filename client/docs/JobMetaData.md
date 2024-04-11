# JobMetaData


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **object** |  | 
**job_id** | **object** |  | 
**target_id** | **object** |  | 
**ligand_id** | **object** |  | 
**folding_method** | **object** |  | 
**docking_method** | **object** |  | 

## Example

```python
from nolabs_microservice.models.job_meta_data import JobMetaData

# TODO update the JSON string below
json = "{}"
# create an instance of JobMetaData from a JSON string
job_meta_data_instance = JobMetaData.from_json(json)
# print the JSON string representation of the object
print JobMetaData.to_json()

# convert the object into a dict
job_meta_data_dict = job_meta_data_instance.to_dict()
# create an instance of JobMetaData from a dict
job_meta_data_form_dict = job_meta_data.from_dict(job_meta_data_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


