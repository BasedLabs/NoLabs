# RunGeneOntologyResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**experiment_id** | **object** |  | 
**experiment_name** | **object** |  | 
**amino_acids** | **object** |  | 

## Example

```python
from nolabs_microservice.models.run_gene_ontology_response import RunGeneOntologyResponse

# TODO update the JSON string below
json = "{}"
# create an instance of RunGeneOntologyResponse from a JSON string
run_gene_ontology_response_instance = RunGeneOntologyResponse.from_json(json)
# print the JSON string representation of the object
print RunGeneOntologyResponse.to_json()

# convert the object into a dict
run_gene_ontology_response_dict = run_gene_ontology_response_instance.to_dict()
# create an instance of RunGeneOntologyResponse from a dict
run_gene_ontology_response_form_dict = run_gene_ontology_response.from_dict(run_gene_ontology_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


