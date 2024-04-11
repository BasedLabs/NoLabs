# nolabs_microservice.FoldingApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**change_experiment_name_api_v1_folding_change_experiment_name_post**](FoldingApi.md#change_experiment_name_api_v1_folding_change_experiment_name_post) | **POST** /api/v1/folding/change-experiment-name | Change Experiment Name
[**create_experiment_api_v1_folding_create_experiment_get**](FoldingApi.md#create_experiment_api_v1_folding_create_experiment_get) | **GET** /api/v1/folding/create-experiment | Create Experiment
[**delete_experiment_api_v1_folding_delete_experiment_delete**](FoldingApi.md#delete_experiment_api_v1_folding_delete_experiment_delete) | **DELETE** /api/v1/folding/delete-experiment | Delete Experiment
[**experiments_api_v1_folding_experiments_get**](FoldingApi.md#experiments_api_v1_folding_experiments_get) | **GET** /api/v1/folding/experiments | Experiments
[**get_experiment_api_v1_folding_get_experiment_get**](FoldingApi.md#get_experiment_api_v1_folding_get_experiment_get) | **GET** /api/v1/folding/get-experiment | Get Experiment
[**inference_api_v1_folding_inference_post**](FoldingApi.md#inference_api_v1_folding_inference_post) | **POST** /api/v1/folding/inference | Inference


# **change_experiment_name_api_v1_folding_change_experiment_name_post**
> object change_experiment_name_api_v1_folding_change_experiment_name_post(change_experiment_name_request)

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
    api_instance = nolabs_microservice.FoldingApi(api_client)
    change_experiment_name_request = nolabs_microservice.ChangeExperimentNameRequest() # ChangeExperimentNameRequest | 

    try:
        # Change Experiment Name
        api_response = api_instance.change_experiment_name_api_v1_folding_change_experiment_name_post(change_experiment_name_request)
        print("The response of FoldingApi->change_experiment_name_api_v1_folding_change_experiment_name_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling FoldingApi->change_experiment_name_api_v1_folding_change_experiment_name_post: %s\n" % e)
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

# **create_experiment_api_v1_folding_create_experiment_get**
> ExperimentMetadataResponse create_experiment_api_v1_folding_create_experiment_get()

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
    api_instance = nolabs_microservice.FoldingApi(api_client)

    try:
        # Create Experiment
        api_response = api_instance.create_experiment_api_v1_folding_create_experiment_get()
        print("The response of FoldingApi->create_experiment_api_v1_folding_create_experiment_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling FoldingApi->create_experiment_api_v1_folding_create_experiment_get: %s\n" % e)
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

# **delete_experiment_api_v1_folding_delete_experiment_delete**
> object delete_experiment_api_v1_folding_delete_experiment_delete(experiment_id)

Delete Experiment

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
    api_instance = nolabs_microservice.FoldingApi(api_client)
    experiment_id = 'experiment_id_example' # str | 

    try:
        # Delete Experiment
        api_response = api_instance.delete_experiment_api_v1_folding_delete_experiment_delete(experiment_id)
        print("The response of FoldingApi->delete_experiment_api_v1_folding_delete_experiment_delete:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling FoldingApi->delete_experiment_api_v1_folding_delete_experiment_delete: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | **str**|  | 

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
**422** | Validation Error |  -  |

[[Back to top]](#) [[Back to API list]](../README.md#documentation-for-api-endpoints) [[Back to Model list]](../README.md#documentation-for-models) [[Back to README]](../README.md)

# **experiments_api_v1_folding_experiments_get**
> object experiments_api_v1_folding_experiments_get()

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
    api_instance = nolabs_microservice.FoldingApi(api_client)

    try:
        # Experiments
        api_response = api_instance.experiments_api_v1_folding_experiments_get()
        print("The response of FoldingApi->experiments_api_v1_folding_experiments_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling FoldingApi->experiments_api_v1_folding_experiments_get: %s\n" % e)
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

# **get_experiment_api_v1_folding_get_experiment_get**
> NolabsApiModelsAminoAcidFoldingGetExperimentResponse get_experiment_api_v1_folding_get_experiment_get(experiment_id)

Get Experiment

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.nolabs_api_models_amino_acid_folding_get_experiment_response import NolabsApiModelsAminoAcidFoldingGetExperimentResponse
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
    api_instance = nolabs_microservice.FoldingApi(api_client)
    experiment_id = 'experiment_id_example' # str | 

    try:
        # Get Experiment
        api_response = api_instance.get_experiment_api_v1_folding_get_experiment_get(experiment_id)
        print("The response of FoldingApi->get_experiment_api_v1_folding_get_experiment_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling FoldingApi->get_experiment_api_v1_folding_get_experiment_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | **str**|  | 

### Return type

[**NolabsApiModelsAminoAcidFoldingGetExperimentResponse**](NolabsApiModelsAminoAcidFoldingGetExperimentResponse.md)

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

# **inference_api_v1_folding_inference_post**
> RunFoldingResponse inference_api_v1_folding_inference_post(experiment_name, experiment_id=experiment_id, amino_acid_sequence=amino_acid_sequence, fastas=fastas)

Inference

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.run_folding_response import RunFoldingResponse
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
    api_instance = nolabs_microservice.FoldingApi(api_client)
    experiment_name = None # object | 
    experiment_id = None # object |  (optional)
    amino_acid_sequence = None # object |  (optional)
    fastas = None # object |  (optional)

    try:
        # Inference
        api_response = api_instance.inference_api_v1_folding_inference_post(experiment_name, experiment_id=experiment_id, amino_acid_sequence=amino_acid_sequence, fastas=fastas)
        print("The response of FoldingApi->inference_api_v1_folding_inference_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling FoldingApi->inference_api_v1_folding_inference_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_name** | [**object**](object.md)|  | 
 **experiment_id** | [**object**](object.md)|  | [optional] 
 **amino_acid_sequence** | [**object**](object.md)|  | [optional] 
 **fastas** | [**object**](object.md)|  | [optional] 

### Return type

[**RunFoldingResponse**](RunFoldingResponse.md)

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

