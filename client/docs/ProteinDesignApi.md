# nolabs_microservice.ProteinDesignApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**change_experiment_name_api_v1_protein_design_change_experiment_name_post**](ProteinDesignApi.md#change_experiment_name_api_v1_protein_design_change_experiment_name_post) | **POST** /api/v1/protein-design/change-experiment-name | Change Experiment Name
[**create_experiment_api_v1_protein_design_create_experiment_get**](ProteinDesignApi.md#create_experiment_api_v1_protein_design_create_experiment_get) | **GET** /api/v1/protein-design/create-experiment | Create Experiment
[**delete_experiment_api_v1_protein_design_experiment_delete**](ProteinDesignApi.md#delete_experiment_api_v1_protein_design_experiment_delete) | **DELETE** /api/v1/protein-design/experiment | Delete Experiment
[**experiments_api_v1_protein_design_experiments_metadata_get**](ProteinDesignApi.md#experiments_api_v1_protein_design_experiments_metadata_get) | **GET** /api/v1/protein-design/experiments-metadata | Experiments
[**get_experiment_api_v1_protein_design_experiment_get**](ProteinDesignApi.md#get_experiment_api_v1_protein_design_experiment_get) | **GET** /api/v1/protein-design/experiment | Get Experiment
[**inference_api_v1_protein_design_inference_post**](ProteinDesignApi.md#inference_api_v1_protein_design_inference_post) | **POST** /api/v1/protein-design/inference | Inference


# **change_experiment_name_api_v1_protein_design_change_experiment_name_post**
> object change_experiment_name_api_v1_protein_design_change_experiment_name_post(change_experiment_name_request)

Change Experiment Name

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.change_experiment_name_request import ChangeExperimentNameRequest
from nolabs_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = nolabs_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with nolabs_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = nolabs_microservice.ProteinDesignApi(api_client)
    change_experiment_name_request = nolabs_microservice.ChangeExperimentNameRequest() # ChangeExperimentNameRequest | 

    try:
        # Change Experiment Name
        api_response = api_instance.change_experiment_name_api_v1_protein_design_change_experiment_name_post(change_experiment_name_request)
        print("The response of ProteinDesignApi->change_experiment_name_api_v1_protein_design_change_experiment_name_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ProteinDesignApi->change_experiment_name_api_v1_protein_design_change_experiment_name_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **change_experiment_name_request** | [**ChangeExperimentNameRequest**](ChangeExperimentNameRequest.md)|  | 

### Return type

**object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: application/json
 - **Accept**: application/json

### HTTP response details

| Status code | Description | Response headers |
|-------------|-------------|------------------|
**200** | Successful Response |  -  |
**422** | Validation Error |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **create_experiment_api_v1_protein_design_create_experiment_get**
> ExperimentMetadataResponse create_experiment_api_v1_protein_design_create_experiment_get()

Create Experiment

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.experiment_metadata_response import ExperimentMetadataResponse
from nolabs_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = nolabs_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with nolabs_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = nolabs_microservice.ProteinDesignApi(api_client)

    try:
        # Create Experiment
        api_response = api_instance.create_experiment_api_v1_protein_design_create_experiment_get()
        print("The response of ProteinDesignApi->create_experiment_api_v1_protein_design_create_experiment_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ProteinDesignApi->create_experiment_api_v1_protein_design_create_experiment_get: %s\n" % e)
```



### Parameters

This endpoint does not need any parameter.

### Return type

[**ExperimentMetadataResponse**](ExperimentMetadataResponse.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details

| Status code | Description | Response headers |
|-------------|-------------|------------------|
**200** | Successful Response |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **delete_experiment_api_v1_protein_design_experiment_delete**
> NolabsApiModelsProteinDesignGetExperimentResponse delete_experiment_api_v1_protein_design_experiment_delete(id)

Delete Experiment

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.nolabs_api_models_protein_design_get_experiment_response import NolabsApiModelsProteinDesignGetExperimentResponse
from nolabs_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = nolabs_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with nolabs_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = nolabs_microservice.ProteinDesignApi(api_client)
    id = None # object | 

    try:
        # Delete Experiment
        api_response = api_instance.delete_experiment_api_v1_protein_design_experiment_delete(id)
        print("The response of ProteinDesignApi->delete_experiment_api_v1_protein_design_experiment_delete:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ProteinDesignApi->delete_experiment_api_v1_protein_design_experiment_delete: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **id** | [**object**](.md)|  | 

### Return type

[**NolabsApiModelsProteinDesignGetExperimentResponse**](NolabsApiModelsProteinDesignGetExperimentResponse.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details

| Status code | Description | Response headers |
|-------------|-------------|------------------|
**200** | Successful Response |  -  |
**422** | Validation Error |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **experiments_api_v1_protein_design_experiments_metadata_get**
> object experiments_api_v1_protein_design_experiments_metadata_get()

Experiments

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = nolabs_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with nolabs_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = nolabs_microservice.ProteinDesignApi(api_client)

    try:
        # Experiments
        api_response = api_instance.experiments_api_v1_protein_design_experiments_metadata_get()
        print("The response of ProteinDesignApi->experiments_api_v1_protein_design_experiments_metadata_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ProteinDesignApi->experiments_api_v1_protein_design_experiments_metadata_get: %s\n" % e)
```



### Parameters

This endpoint does not need any parameter.

### Return type

**object**

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details

| Status code | Description | Response headers |
|-------------|-------------|------------------|
**200** | Successful Response |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **get_experiment_api_v1_protein_design_experiment_get**
> NolabsApiModelsProteinDesignGetExperimentResponse get_experiment_api_v1_protein_design_experiment_get(id)

Get Experiment

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.nolabs_api_models_protein_design_get_experiment_response import NolabsApiModelsProteinDesignGetExperimentResponse
from nolabs_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = nolabs_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with nolabs_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = nolabs_microservice.ProteinDesignApi(api_client)
    id = None # object | 

    try:
        # Get Experiment
        api_response = api_instance.get_experiment_api_v1_protein_design_experiment_get(id)
        print("The response of ProteinDesignApi->get_experiment_api_v1_protein_design_experiment_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ProteinDesignApi->get_experiment_api_v1_protein_design_experiment_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **id** | [**object**](.md)|  | 

### Return type

[**NolabsApiModelsProteinDesignGetExperimentResponse**](NolabsApiModelsProteinDesignGetExperimentResponse.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: Not defined
 - **Accept**: application/json

### HTTP response details

| Status code | Description | Response headers |
|-------------|-------------|------------------|
**200** | Successful Response |  -  |
**422** | Validation Error |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **inference_api_v1_protein_design_inference_post**
> RunProteinDesignResponse inference_api_v1_protein_design_inference_post(experiment_name, pdb_file, experiment_id=experiment_id, contig=contig, number_of_designs=number_of_designs, timesteps=timesteps, hotspots=hotspots)

Inference

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.run_protein_design_response import RunProteinDesignResponse
from nolabs_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = nolabs_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with nolabs_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = nolabs_microservice.ProteinDesignApi(api_client)
    experiment_name = None # object | 
    pdb_file = None # object | 
    experiment_id = None # object |  (optional)
    contig = None # object |  (optional)
    number_of_designs = None # object |  (optional)
    timesteps = None # object |  (optional)
    hotspots = None # object |  (optional)

    try:
        # Inference
        api_response = api_instance.inference_api_v1_protein_design_inference_post(experiment_name, pdb_file, experiment_id=experiment_id, contig=contig, number_of_designs=number_of_designs, timesteps=timesteps, hotspots=hotspots)
        print("The response of ProteinDesignApi->inference_api_v1_protein_design_inference_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ProteinDesignApi->inference_api_v1_protein_design_inference_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_name** | [**object**](object.md)|  | 
 **pdb_file** | [**object**](object.md)|  | 
 **experiment_id** | [**object**](object.md)|  | [optional] 
 **contig** | [**object**](object.md)|  | [optional] 
 **number_of_designs** | [**object**](object.md)|  | [optional] 
 **timesteps** | [**object**](object.md)|  | [optional] 
 **hotspots** | [**object**](object.md)|  | [optional] 

### Return type

[**RunProteinDesignResponse**](RunProteinDesignResponse.md)

### Authorization

No authorization required

### HTTP request headers

 - **Content-Type**: multipart/form-data
 - **Accept**: application/json

### HTTP response details

| Status code | Description | Response headers |
|-------------|-------------|------------------|
**200** | Successful Response |  -  |
**422** | Validation Error |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

