# rcsb_pdb_query_microservice.DefaultApi

All URIs are relative to *http://localhost*

Method | HTTP request | Description
------------- | ------------- | -------------
[**get_running_jobs_jobs_running_get**](DefaultApi.md#get_running_jobs_jobs_running_get) | **GET** /jobs/running | Get Running Jobs
[**is_job_running_job_job_id_is_running_get**](DefaultApi.md#is_job_running_job_job_id_is_running_get) | **GET** /job/{job_id}/is-running | Is Job Running
[**predict_fetch_fastas_by_ids_post**](DefaultApi.md#predict_fetch_fastas_by_ids_post) | **POST** /fetch-fastas-by-ids | Predict


# **get_running_jobs_jobs_running_get**
> object get_running_jobs_jobs_running_get()

Get Running Jobs

### Example


```python
import rcsb_pdb_query_microservice
from rcsb_pdb_query_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = rcsb_pdb_query_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with rcsb_pdb_query_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = rcsb_pdb_query_microservice.DefaultApi(api_client)

    try:
        # Get Running Jobs
        api_response = api_instance.get_running_jobs_jobs_running_get()
        print("The response of DefaultApi->get_running_jobs_jobs_running_get:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->get_running_jobs_jobs_running_get: %s\n" % e)
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

# **is_job_running_job_job_id_is_running_get**
> IsJobRunningResponse is_job_running_job_job_id_is_running_get(job_id)

Is Job Running

### Example


```python
import rcsb_pdb_query_microservice
from rcsb_pdb_query_microservice.models.is_job_running_response import IsJobRunningResponse
from rcsb_pdb_query_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = rcsb_pdb_query_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with rcsb_pdb_query_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = rcsb_pdb_query_microservice.DefaultApi(api_client)
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

# **predict_fetch_fastas_by_ids_post**
> GetFastaFilesByIdsResponse predict_fetch_fastas_by_ids_post(get_fasta_files_by_ids_request)

Predict

### Example


```python
import rcsb_pdb_query_microservice
from rcsb_pdb_query_microservice.models.get_fasta_files_by_ids_request import GetFastaFilesByIdsRequest
from rcsb_pdb_query_microservice.models.get_fasta_files_by_ids_response import GetFastaFilesByIdsResponse
from rcsb_pdb_query_microservice.rest import ApiException
from pprint import pprint

# Defining the host is optional and defaults to http://localhost
# See configuration.py for a list of all supported configuration parameters.
configuration = rcsb_pdb_query_microservice.Configuration(
    host = "http://localhost"
)


# Enter a context with an instance of the API client
with rcsb_pdb_query_microservice.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = rcsb_pdb_query_microservice.DefaultApi(api_client)
    get_fasta_files_by_ids_request = rcsb_pdb_query_microservice.GetFastaFilesByIdsRequest() # GetFastaFilesByIdsRequest | 

    try:
        # Predict
        api_response = api_instance.predict_fetch_fastas_by_ids_post(get_fasta_files_by_ids_request)
        print("The response of DefaultApi->predict_fetch_fastas_by_ids_post:\n")
        pprint(api_response)
    except Exception as e:
        print("Exception when calling DefaultApi->predict_fetch_fastas_by_ids_post: %s\n" % e)
```



### Parameters


Name | Type | Description  | Notes
------------- | ------------- | ------------- | -------------
 **get_fasta_files_by_ids_request** | [**GetFastaFilesByIdsRequest**](GetFastaFilesByIdsRequest.md)|  | 

### Return type

[**GetFastaFilesByIdsResponse**](GetFastaFilesByIdsResponse.md)

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

