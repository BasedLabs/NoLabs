# DrugIndicationRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**filters** | **object** |  | [optional] 
**order_by** | **str** |  | [optional] 
**limit** | **int** |  | [optional] [default to 20]
**job_id** | [**JobId**](JobId.md) |  | [optional] 

## Example

```python
from chembl_query_microservice.models.drug_indication_request import DrugIndicationRequest

# TODO update the JSON string below
json = "{}"
# create an instance of DrugIndicationRequest from a JSON string
drug_indication_request_instance = DrugIndicationRequest.from_json(json)
# print the JSON string representation of the object
print DrugIndicationRequest.to_json()

# convert the object into a dict
drug_indication_request_dict = drug_indication_request_instance.to_dict()
# create an instance of DrugIndicationRequest from a dict
drug_indication_request_form_dict = drug_indication_request.from_dict(drug_indication_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


