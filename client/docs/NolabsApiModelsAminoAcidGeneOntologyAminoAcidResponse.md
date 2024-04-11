# NolabsApiModelsAminoAcidGeneOntologyAminoAcidResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**sequence** | **object** |  | 
**name** | **object** |  | 
**go** | [**Dict[str, RunGeneOntologyResponseDataNode]**](RunGeneOntologyResponseDataNode.md) |  | 

## Example

```python
from nolabs_microservice.models.nolabs_api_models_amino_acid_gene_ontology_amino_acid_response import NolabsApiModelsAminoAcidGeneOntologyAminoAcidResponse

# TODO update the JSON string below
json = "{}"
# create an instance of NolabsApiModelsAminoAcidGeneOntologyAminoAcidResponse from a JSON string
nolabs_api_models_amino_acid_gene_ontology_amino_acid_response_instance = NolabsApiModelsAminoAcidGeneOntologyAminoAcidResponse.from_json(json)
# print the JSON string representation of the object
print NolabsApiModelsAminoAcidGeneOntologyAminoAcidResponse.to_json()

# convert the object into a dict
nolabs_api_models_amino_acid_gene_ontology_amino_acid_response_dict = nolabs_api_models_amino_acid_gene_ontology_amino_acid_response_instance.to_dict()
# create an instance of NolabsApiModelsAminoAcidGeneOntologyAminoAcidResponse from a dict
nolabs_api_models_amino_acid_gene_ontology_amino_acid_response_form_dict = nolabs_api_models_amino_acid_gene_ontology_amino_acid_response.from_dict(nolabs_api_models_amino_acid_gene_ontology_amino_acid_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


