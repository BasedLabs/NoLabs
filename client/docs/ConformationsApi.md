# nolabs_microservice.ConformationsApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**change_experiment_name_api_v1_conformations_change_experiment_name_post**](ConformationsApi.md#change_experiment_name_api_v1_conformations_change_experiment_name_post) | **POST** /api/v1/conformations/change-experiment-name | Change Experiment Name
[**create_experiment_api_v1_conformations_create_experiment_get**](ConformationsApi.md#create_experiment_api_v1_conformations_create_experiment_get) | **GET** /api/v1/conformations/create-experiment | Create Experiment
[**delete_experiment_api_v1_conformations_delete_experiment_delete**](ConformationsApi.md#delete_experiment_api_v1_conformations_delete_experiment_delete) | **DELETE** /api/v1/conformations/delete-experiment | Delete Experiment
[**experiments_api_v1_conformations_experiments_get**](ConformationsApi.md#experiments_api_v1_conformations_experiments_get) | **GET** /api/v1/conformations/experiments | Experiments
[**get_experiment_api_v1_conformations_get_experiment_get**](ConformationsApi.md#get_experiment_api_v1_conformations_get_experiment_get) | **GET** /api/v1/conformations/get-experiment | Get Experiment
[**inference_api_v1_conformations_inference_post**](ConformationsApi.md#inference_api_v1_conformations_inference_post) | **POST** /api/v1/conformations/inference | Inference


# **change_experiment_name_api_v1_conformations_change_experiment_name_post**
> object change_experiment_name_api_v1_conformations_change_experiment_name_post(change_experiment_name_request)

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
    api_instance = nolabs_microservice.ConformationsApi(api_client)
    change_experiment_name_request = nolabs_microservice.ChangeExperimentNameRequest() # ChangeExperimentNameRequest | 

    try:
        # Change Experiment Name
        api_response = api_instance.change_experiment_name_api_v1_conformations_change_experiment_name_post(change_experiment_name_request)
        print("The response of ConformationsApi->change_experiment_name_api_v1_conformations_change_experiment_name_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ConformationsApi->change_experiment_name_api_v1_conformations_change_experiment_name_post: %s\n" % e)
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

# **create_experiment_api_v1_conformations_create_experiment_get**
> ExperimentMetadataResponse create_experiment_api_v1_conformations_create_experiment_get()

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
    api_instance = nolabs_microservice.ConformationsApi(api_client)

    try:
        # Create Experiment
        api_response = api_instance.create_experiment_api_v1_conformations_create_experiment_get()
        print("The response of ConformationsApi->create_experiment_api_v1_conformations_create_experiment_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ConformationsApi->create_experiment_api_v1_conformations_create_experiment_get: %s\n" % e)
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

# **delete_experiment_api_v1_conformations_delete_experiment_delete**
> object delete_experiment_api_v1_conformations_delete_experiment_delete(experiment_id)

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
    api_instance = nolabs_microservice.ConformationsApi(api_client)
    experiment_id = 'experiment_id_example' # str | 

    try:
        # Delete Experiment
        api_response = api_instance.delete_experiment_api_v1_conformations_delete_experiment_delete(experiment_id)
        print("The response of ConformationsApi->delete_experiment_api_v1_conformations_delete_experiment_delete:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ConformationsApi->delete_experiment_api_v1_conformations_delete_experiment_delete: %s\n" % e)
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

# **experiments_api_v1_conformations_experiments_get**
> List[ExperimentMetadataResponse] experiments_api_v1_conformations_experiments_get()

Experiments

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
    api_instance = nolabs_microservice.ConformationsApi(api_client)

    try:
        # Experiments
        api_response = api_instance.experiments_api_v1_conformations_experiments_get()
        print("The response of ConformationsApi->experiments_api_v1_conformations_experiments_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ConformationsApi->experiments_api_v1_conformations_experiments_get: %s\n" % e)
```



### Parameters

This endpoint does not need any parameter.

### Return type

[**List[ExperimentMetadataResponse]**](ExperimentMetadataResponse.md)

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

# **get_experiment_api_v1_conformations_get_experiment_get**
> NolabsApiModelsConformationsGetExperimentResponse get_experiment_api_v1_conformations_get_experiment_get(experiment_id)

Get Experiment

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.nolabs_api_models_conformations_get_experiment_response import NolabsApiModelsConformationsGetExperimentResponse
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
    api_instance = nolabs_microservice.ConformationsApi(api_client)
    experiment_id = 'experiment_id_example' # str | 

    try:
        # Get Experiment
        api_response = api_instance.get_experiment_api_v1_conformations_get_experiment_get(experiment_id)
        print("The response of ConformationsApi->get_experiment_api_v1_conformations_get_experiment_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ConformationsApi->get_experiment_api_v1_conformations_get_experiment_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **experiment_id** | **str**|  | 

### Return type

[**NolabsApiModelsConformationsGetExperimentResponse**](NolabsApiModelsConformationsGetExperimentResponse.md)

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

# **inference_api_v1_conformations_inference_post**
> RunSimulationsResponse inference_api_v1_conformations_inference_post(pdb_file, experiment_name, experiment_id=experiment_id, total_frames=total_frames, temperature_k=temperature_k, take_frame_every=take_frame_every, step_size=step_size, replace_non_standard_residues=replace_non_standard_residues, add_missing_atoms=add_missing_atoms, add_missing_hydrogens=add_missing_hydrogens, friction_coeff=friction_coeff, ignore_missing_atoms=ignore_missing_atoms, integrator=integrator)

Inference

### Example


```python
import time
import os
import nolabs_microservice
from nolabs_microservice.models.run_simulations_response import RunSimulationsResponse
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
    api_instance = nolabs_microservice.ConformationsApi(api_client)
    pdb_file = None # object | 
    experiment_name = None # object | 
    experiment_id = None # object |  (optional)
    total_frames = None # object |  (optional)
    temperature_k = None # object |  (optional)
    take_frame_every = None # object |  (optional)
    step_size = None # object |  (optional)
    replace_non_standard_residues = None # object |  (optional)
    add_missing_atoms = None # object |  (optional)
    add_missing_hydrogens = None # object |  (optional)
    friction_coeff = None # object |  (optional)
    ignore_missing_atoms = None # object |  (optional)
    integrator = nolabs_microservice.IntegratorsRequest() # IntegratorsRequest |  (optional)

    try:
        # Inference
        api_response = api_instance.inference_api_v1_conformations_inference_post(pdb_file, experiment_name, experiment_id=experiment_id, total_frames=total_frames, temperature_k=temperature_k, take_frame_every=take_frame_every, step_size=step_size, replace_non_standard_residues=replace_non_standard_residues, add_missing_atoms=add_missing_atoms, add_missing_hydrogens=add_missing_hydrogens, friction_coeff=friction_coeff, ignore_missing_atoms=ignore_missing_atoms, integrator=integrator)
        print("The response of ConformationsApi->inference_api_v1_conformations_inference_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling ConformationsApi->inference_api_v1_conformations_inference_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **pdb_file** | [**object**](object.md)|  | 
 **experiment_name** | [**object**](object.md)|  | 
 **experiment_id** | [**object**](object.md)|  | [optional] 
 **total_frames** | [**object**](object.md)|  | [optional] 
 **temperature_k** | [**object**](object.md)|  | [optional] 
 **take_frame_every** | [**object**](object.md)|  | [optional] 
 **step_size** | [**object**](object.md)|  | [optional] 
 **replace_non_standard_residues** | [**object**](object.md)|  | [optional] 
 **add_missing_atoms** | [**object**](object.md)|  | [optional] 
 **add_missing_hydrogens** | [**object**](object.md)|  | [optional] 
 **friction_coeff** | [**object**](object.md)|  | [optional] 
 **ignore_missing_atoms** | [**object**](object.md)|  | [optional] 
 **integrator** | [**IntegratorsRequest**](IntegratorsRequest.md)|  | [optional] 

### Return type

[**RunSimulationsResponse**](RunSimulationsResponse.md)

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

