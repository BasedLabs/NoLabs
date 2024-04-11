# RunLocalisationResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **object** |  | 
**experiment_name** | **object** |  | 
**amino_acids** | **object** |  | 

## Example

```python
from nolabs_microservice.models.run_localisation_response import RunLocalisationResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunLocalisationResponse from a JSON string
run_localisation_response_instance = RunLocalisationResponse.from_json(json)
# print the JSON string representation of the object
print RunLocalisationResponse.to_json()

# convert the object into a dict
run_localisation_response_dict = run_localisation_response_instance.to_dict()
# create an instance of RunLocalisationResponse from a dict
run_localisation_response_form_dict = run_localisation_response.from_dict(run_localisation_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


