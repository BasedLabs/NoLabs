# RunGeneOntologyPredictionRequest


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**amino_acid_sequence** | **str** |  | 

## Example

```python
from gene_ontology_microservice.models.run_gene_ontology_prediction_request import RunGeneOntologyPredictionRequest

# TODO update the JSON string below
json = "{}"
# create an instance of RunGeneOntologyPredictionRequest from a JSON string
run_gene_ontology_prediction_request_instance = RunGeneOntologyPredictionRequest.from_json(json)
# print the JSON string representation of the object
print RunGeneOntologyPredictionRequest.to_json()

# convert the object into a dict
run_gene_ontology_prediction_request_dict = run_gene_ontology_prediction_request_instance.to_dict()
# create an instance of RunGeneOntologyPredictionRequest from a dict
run_gene_ontology_prediction_request_form_dict = run_gene_ontology_prediction_request.from_dict(run_gene_ontology_prediction_request_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


