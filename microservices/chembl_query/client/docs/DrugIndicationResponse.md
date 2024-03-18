# DrugIndicationResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**drugs** | [**List[Molecule]**](Molecule.md) |  | 
**total_count** | **int** |  | 
**page** | **int** |  | 
**pages** | **int** |  | 

## Example

```python
from chembl_query_microservice.models.drug_indication_response import DrugIndicationResponse

# TODO update the JSON string below
json = "{}"
# create an instance of DrugIndicationResponse from a JSON string
drug_indication_response_instance = DrugIndicationResponse.from_json(json)
# print the JSON string representation of the object
print DrugIndicationResponse.to_json()

# convert the object into a dict
drug_indication_response_dict = drug_indication_response_instance.to_dict()
# create an instance of DrugIndicationResponse from a dict
drug_indication_response_form_dict = drug_indication_response.from_dict(drug_indication_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


