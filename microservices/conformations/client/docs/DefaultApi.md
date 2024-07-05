# conformations_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**gen_gro_top_endpoint_gen_gro_top_post**](DefaultApi.md#gen_gro_top_endpoint_gen_gro_top_post) | **POST** /gen-gro-top | Gen Gro Top Endpoint
[**is_job_running_job_job_id_is_running_get**](DefaultApi.md#is_job_running_job_job_id_is_running_get) | **GET** /job/{job_id}/is-running | Is Job Running
[**run_gromacs_simulations_endpoint_run_gromacs_simulations_post**](DefaultApi.md#run_gromacs_simulations_endpoint_run_gromacs_simulations_post) | **POST** /run-gromacs-simulations | Run Gromacs Simulations Endpoint
[**run_pdb_fixer_endpoint_run_pdb_fixer_post**](DefaultApi.md#run_pdb_fixer_endpoint_run_pdb_fixer_post) | **POST** /run-pdb-fixer | Run Pdb Fixer Endpoint
[**run_pdb_simulations_endpoint_run_pdb_simulations_post**](DefaultApi.md#run_pdb_simulations_endpoint_run_pdb_simulations_post) | **POST** /run-pdb-simulations | Run Pdb Simulations Endpoint


# **gen_gro_top_endpoint_gen_gro_top_post**
> GenGroTopResponse gen_gro_top_endpoint_gen_gro_top_post(gen_gro_top_request)

Gen Gro Top Endpoint

### Example


```python
import time
import os
import conformations_microservice
from conformations_microservice.models.gen_gro_top_request import GenGroTopRequest
from conformations_microservice.models.gen_gro_top_response import GenGroTopResponse
from conformations_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = conformations_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with conformations_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = conformations_microservice.DefaultApi(api_client)
    gen_gro_top_request = conformations_microservice.GenGroTopRequest() # GenGroTopRequest | 

    try:
        # Gen Gro Top Endpoint
        api_response = api_instance.gen_gro_top_endpoint_gen_gro_top_post(gen_gro_top_request)
        print("The response of DefaultApi->gen_gro_top_endpoint_gen_gro_top_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->gen_gro_top_endpoint_gen_gro_top_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **gen_gro_top_request** | [**GenGroTopRequest**](GenGroTopRequest.md)|  | 

### Return type

[**GenGroTopResponse**](GenGroTopResponse.md)

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

# **is_job_running_job_job_id_is_running_get**
> IsJobRunningResponse is_job_running_job_job_id_is_running_get(job_id)

Is Job Running

### Example


```python
import time
import os
import conformations_microservice
from conformations_microservice.models.is_job_running_response import IsJobRunningResponse
from conformations_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = conformations_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with conformations_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = conformations_microservice.DefaultApi(api_client)
    job_id = 'job_id_example' # str | 

    try:
        # Is Job Running
        api_response = api_instance.is_job_running_job_job_id_is_running_get(job_id)
        print("The response of DefaultApi->is_job_running_job_job_id_is_running_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->is_job_running_job_job_id_is_running_get: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | **str**|  | 

### Return type

[**IsJobRunningResponse**](IsJobRunningResponse.md)

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

# **run_gromacs_simulations_endpoint_run_gromacs_simulations_post**
> RunSimulationsResponse run_gromacs_simulations_endpoint_run_gromacs_simulations_post(run_gromacs_simulations_request)

Run Gromacs Simulations Endpoint

### Example


```python
import time
import os
import conformations_microservice
from conformations_microservice.models.run_gromacs_simulations_request import RunGromacsSimulationsRequest
from conformations_microservice.models.run_simulations_response import RunSimulationsResponse
from conformations_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = conformations_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with conformations_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = conformations_microservice.DefaultApi(api_client)
    run_gromacs_simulations_request = conformations_microservice.RunGromacsSimulationsRequest() # RunGromacsSimulationsRequest | 

    try:
        # Run Gromacs Simulations Endpoint
        api_response = api_instance.run_gromacs_simulations_endpoint_run_gromacs_simulations_post(run_gromacs_simulations_request)
        print("The response of DefaultApi->run_gromacs_simulations_endpoint_run_gromacs_simulations_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->run_gromacs_simulations_endpoint_run_gromacs_simulations_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **run_gromacs_simulations_request** | [**RunGromacsSimulationsRequest**](RunGromacsSimulationsRequest.md)|  | 

### Return type

[**RunSimulationsResponse**](RunSimulationsResponse.md)

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

# **run_pdb_fixer_endpoint_run_pdb_fixer_post**
> RunPdbFixerResponse run_pdb_fixer_endpoint_run_pdb_fixer_post(run_pdb_fixer_request)

Run Pdb Fixer Endpoint

### Example


```python
import time
import os
import conformations_microservice
from conformations_microservice.models.run_pdb_fixer_request import RunPdbFixerRequest
from conformations_microservice.models.run_pdb_fixer_response import RunPdbFixerResponse
from conformations_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = conformations_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with conformations_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = conformations_microservice.DefaultApi(api_client)
    run_pdb_fixer_request = conformations_microservice.RunPdbFixerRequest() # RunPdbFixerRequest | 

    try:
        # Run Pdb Fixer Endpoint
        api_response = api_instance.run_pdb_fixer_endpoint_run_pdb_fixer_post(run_pdb_fixer_request)
        print("The response of DefaultApi->run_pdb_fixer_endpoint_run_pdb_fixer_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->run_pdb_fixer_endpoint_run_pdb_fixer_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **run_pdb_fixer_request** | [**RunPdbFixerRequest**](RunPdbFixerRequest.md)|  | 

### Return type

[**RunPdbFixerResponse**](RunPdbFixerResponse.md)

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

# **run_pdb_simulations_endpoint_run_pdb_simulations_post**
> RunSimulationsResponse run_pdb_simulations_endpoint_run_pdb_simulations_post(run_pdb_simulations_request)

Run Pdb Simulations Endpoint

### Example


```python
import time
import os
import conformations_microservice
from conformations_microservice.models.run_pdb_simulations_request import RunPdbSimulationsRequest
from conformations_microservice.models.run_simulations_response import RunSimulationsResponse
from conformations_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = conformations_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with conformations_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = conformations_microservice.DefaultApi(api_client)
    run_pdb_simulations_request = conformations_microservice.RunPdbSimulationsRequest() # RunPdbSimulationsRequest | 

    try:
        # Run Pdb Simulations Endpoint
        api_response = api_instance.run_pdb_simulations_endpoint_run_pdb_simulations_post(run_pdb_simulations_request)
        print("The response of DefaultApi->run_pdb_simulations_endpoint_run_pdb_simulations_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->run_pdb_simulations_endpoint_run_pdb_simulations_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **run_pdb_simulations_request** | [**RunPdbSimulationsRequest**](RunPdbSimulationsRequest.md)|  | 

### Return type

[**RunSimulationsResponse**](RunSimulationsResponse.md)

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

