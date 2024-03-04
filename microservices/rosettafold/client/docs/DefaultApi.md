# rosettafold_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**is_job_running_job_job_id_is_running_get**](DefaultApi.md#is_job_running_job_job_id_is_running_get) | **GET** /job/{job_id}/is-running | Is Job Running
[**livez_livez_get**](DefaultApi.md#livez_livez_get) | **GET** /livez | Livez
[**run_folding_run_folding_post**](DefaultApi.md#run_folding_run_folding_post) | **POST** /run-folding | Run Folding


# **is_job_running_job_job_id_is_running_get**
> IsJobRunningResponse is_job_running_job_job_id_is_running_get(job_id)

Is Job Running

### Example


```python
import rosettafold_microservice
from rosettafold_microservice.models.is_job_running_response import IsJobRunningResponse
from rosettafold_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = rosettafold_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with rosettafold_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = rosettafold_microservice.DefaultApi(api_client)
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

# **livez_livez_get**
> object livez_livez_get()

Livez

### Example


```python
import rosettafold_microservice
from rosettafold_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = rosettafold_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with rosettafold_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = rosettafold_microservice.DefaultApi(api_client)

    try:
        # Livez
        api_response = api_instance.livez_livez_get()
        print("The response of DefaultApi->livez_livez_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->livez_livez_get: %s\n" % e)
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

# **run_folding_run_folding_post**
> RunRosettaFoldResponse run_folding_run_folding_post(job_id=job_id, fasta=fasta, a3m=a3m)

Run Folding

Run folding on a given amino-acid sequence :param job_id str | None: Used for job tracking :param fasta Annotated[bytes, File()]: Fasta file with amino acid sequence. You must specify either fasta or a3m :param a3m Annotated[bytes, File()]: MSA file. You must specify either fasta or a3m

### Example


```python
import rosettafold_microservice
from rosettafold_microservice.models.run_rosetta_fold_response import RunRosettaFoldResponse
from rosettafold_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = rosettafold_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with rosettafold_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = rosettafold_microservice.DefaultApi(api_client)
    job_id = rosettafold_microservice.JobId() # JobId |  (optional)
    fasta = None # bytearray |  (optional)
    a3m = None # bytearray |  (optional)

    try:
        # Run Folding
        api_response = api_instance.run_folding_run_folding_post(job_id=job_id, fasta=fasta, a3m=a3m)
        print("The response of DefaultApi->run_folding_run_folding_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->run_folding_run_folding_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **job_id** | [**JobId**](.md)|  | [optional] 
 **fasta** | **bytearray**|  | [optional] 
 **a3m** | **bytearray**|  | [optional] 

### Return type

[**RunRosettaFoldResponse**](RunRosettaFoldResponse.md)

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

