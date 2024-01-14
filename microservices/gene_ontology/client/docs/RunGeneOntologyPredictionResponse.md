# RunGeneOntologyPredictionResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**go_confidence** | [**List[GoConfidenceResponse]**](GoConfidenceResponse.md) |  | 
**errors** | **List[str]** |  | 

## Example

```python
from gene_ontology_microservice.models.run_gene_ontology_prediction_response import RunGeneOntologyPredictionResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunGeneOntologyPredictionResponse from a JSON string
run_gene_ontology_prediction_response_instance = RunGeneOntologyPredictionResponse.from_json(json)
# print the JSON string representation of the object
print RunGeneOntologyPredictionResponse.to_json()

# convert the object into a dict
run_gene_ontology_prediction_response_dict = run_gene_ontology_prediction_response_instance.to_dict()
# create an instance of RunGeneOntologyPredictionResponse from a dict
run_gene_ontology_prediction_response_form_dict = run_gene_ontology_prediction_response.from_dict(run_gene_ontology_prediction_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


