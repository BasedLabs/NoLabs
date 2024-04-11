# NolabsApiModelsAminoAcidLocalisationAminoAcidResponse


## Properties

Name | Type | Description | Notes
------------ | ------------- | ------------- | -------------
**sequence** | **object** |  | 
**name** | **object** |  | 
**cytosolic_proteins** | **object** |  | 
**mitochondrial_proteins** | **object** |  | 
**nuclear_proteins** | **object** |  | 
**other_proteins** | **object** |  | 
**extracellular_secreted_proteins** | **object** |  | 

## Example

```python
from nolabs_microservice.models.nolabs_api_models_amino_acid_localisation_amino_acid_response import NolabsApiModelsAminoAcidLocalisationAminoAcidResponse

# TODO update the JSON string below
json = "{}"
# create an instance of NolabsApiModelsAminoAcidLocalisationAminoAcidResponse from a JSON string
nolabs_api_models_amino_acid_localisation_amino_acid_response_instance = NolabsApiModelsAminoAcidLocalisationAminoAcidResponse.from_json(json)
# print the JSON string representation of the object
print NolabsApiModelsAminoAcidLocalisationAminoAcidResponse.to_json()

# convert the object into a dict
nolabs_api_models_amino_acid_localisation_amino_acid_response_dict = nolabs_api_models_amino_acid_localisation_amino_acid_response_instance.to_dict()
# create an instance of NolabsApiModelsAminoAcidLocalisationAminoAcidResponse from a dict
nolabs_api_models_amino_acid_localisation_amino_acid_response_form_dict = nolabs_api_models_amino_acid_localisation_amino_acid_response.from_dict(nolabs_api_models_amino_acid_localisation_amino_acid_response_dict)
```
[[Back to Model list]](../README.md#documentation-for-models) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to README]](../README.md)


